//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>
#include <stdlib.h>

#include <cilk/cilk.h>
#include "bigfloat.h"
#include "fft.h"

#include <math.h>

// Multiplies by a uint32_t inline.
void bigfloat_mului(bigfloat_t a, uint32_t b) {
  // Make sure there is enough mem for the result (because we modify the coef
  // pointer and the realloc goes BOOM). If there is not enough mem, the result
  // will be truncated or just be dead wrong.
  //bigfloat_realloc(a, a->len+1);
  a->coef[a->len-1] = 0;

  uint32_t carry = 0;
  for(size_t i=0;i<a->len;i++) {
    uint64_t product = (uint64_t)a->coef[i] * (uint64_t)b + carry;

    uint32_t lower = product & 0xffffffff;
    carry = product >> 32;

    a->coef[i] = lower;
  }
}

// Multiplies by a uint32_t with target.
void bigfloat_mulu(bigfloat_t target, const bigfloat_t a, uint32_t b) {
  // Make sure there is enough mem for the result (because we modify the coef
  // pointer and the realloc goes BOOM). If there is not enough mem, the result
  // will be truncated or just be dead wrong.
  //bigfloat_realloc(a, a->len+1);

  uint32_t carry = 0;
  for(size_t i=0;i<target->len;i++) {
    uint64_t product = carry;

    if(i < a->len) {
      product += (uint64_t)a->coef[i] * (uint64_t)b;
    }

    uint32_t lower = product & 0xffffffff;
    carry = product >> 32;

    target->coef[i] = lower;
  }
}

// Multiplies C = A*B using FFT where C = CL integers in CT, etc.
// CT, CL = target ptr and length
// AT, AL = source 1 ptr and length
// BT, BL = source 2 ptr and length
int static inline _fft_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
  int bits_per_point = 20;

  // Experimentally and guessed
  if(CL > 2000)       bits_per_point = 19;
  if(CL > 9000)       bits_per_point = 18;
  if(CL > 25000)      bits_per_point = 17;
  if(CL > 80000)      bits_per_point = 16;
  if(CL > 280000)     bits_per_point = 15;
  if(CL > 1000000)    bits_per_point = 14;
  if(CL > 4000000)    bits_per_point = 13;
  if(CL > 13000000)   bits_per_point = 12;
  if(CL > 50000000)   bits_per_point = 11;
  if(CL > 200000000)  bits_per_point = 10;
  if(CL > 800000000)  bits_per_point = 9;

  //  Determine minimum FFT size.
  size_t length = fft_length(32.0/bits_per_point*(double)CL);

  // If the arguments are the same, we skip one conversion.
  char needB = (AT != BT || AL != BL);

  //  Allocate FFT arrays
  complex double *Ta = (complex double*)_mm_malloc(length * sizeof(complex double), 32);
  complex double *Tb = NULL;

  if(needB) {
    Tb = (complex double*)_mm_malloc(length * sizeof(complex double), 32);
  }

  //  Convert Numbers to FFT
  int_to_fft(Ta, length, AT, AL, bits_per_point);

  if(needB) {
    int_to_fft(Tb, length, BT, BL, bits_per_point);
  }

  if(needB) {
    cilk_spawn fft_forward(Ta, length);
    fft_forward(Tb, length);
    cilk_sync;
  } else {
    fft_forward(Ta, length);
  }

  // Pointwise multiply
  if(needB) {
    fft_pointwise(Ta, Tb, length);
  } else {
    fft_pointwise(Ta, Ta, length);
  }

  fft_inverse(Ta, length);

  // Convert including carryout
  fft_to_int(Ta, length, CT, CL, bits_per_point);

  // Free FFT arrays
  if(needB) {
    _mm_free(Tb);
  }
  _mm_free(Ta);

  return bits_per_point;
}

void static inline _basecase_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
  char inplaceOverride = (CT == AT || CT == BT);

  if(inplaceOverride) {
    uint32_t* PT = malloc(CL*sizeof(uint32_t));
    _basecase_mul(PT, CL, AT, AL, BT, BL);
    memcpy(CT, PT, CL*sizeof(uint32_t));
    free(PT);
    return;
  }

  // Zero target
  memset(CT, 0, sizeof(uint32_t)*CL);

  for(size_t i=0;i<AL;i++) {
    uint64_t carry = 0;
    uint64_t value;
    for(size_t j=0;j<BL;j++) {
      value = CT[i+j];
      value += (uint64_t) AT[i] * (uint64_t) BT[j];
      value += carry;

      carry  = value >> 32;
      value &= 0xffffffff;
      CT[i+j] = value;
    }
    for(size_t j=i+BL;j<CL;j++) {
      if(carry == 0) {
        break;
      }
      value = CT[j];
      value += carry;

      carry  = value >> 32;
      value &= 0xffffffff;
      CT[j] = value;
    }
  }
}

#define CHECKPRIME (0x9f063b05)

uint32_t checksum(uint32_t *AT, size_t AL) {
  uint64_t base = 1;
  uint64_t prime = CHECKPRIME;
  uint64_t sum = 0;
  for(size_t i=0;i<AL;i++) {
    uint64_t v = AT[i];
    uint64_t product = (v*base) % prime;
    base = (base << 32)%prime;
    sum = sum + product;
    sum %= prime;
  }

  return (uint32_t)sum;
}

uint32_t checksum_mul(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p  = (av*bv)%CHECKPRIME;
  return (uint32_t)p;
}

#define BASECASE_THRESH 1500

void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
    //  Multiplication

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.
    #ifdef DEBUG
    printf("mul a*b = t\n");
    bigfloat_print("a", a);
    bigfloat_print("b", b);
    #endif

    //  Either operand is zero.
    if (bigfloat_iszero(a) || bigfloat_iszero(b)) {
      bigfloat_zero(target);
      return;
    }

    if (p == 0) {
        //  Default value. No trunction.
        p = a->len + b->len;
    } else {
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    //  Collect operands.
    int64_t Aexp = a->exp;
    int64_t Bexp = b->exp;
    size_t AL = a->len;
    size_t BL = b->len;
    uint32_t *AT = a->coef;
    uint32_t *BT = b->coef;

    //  Perform precision truncation.
    if (AL > p) {
        size_t chop = AL - p;
        AL = p;
        Aexp += chop;
        AT += chop;
    }
    if (BL > p) {
        size_t chop = BL - p;
        BL = p;
        Bexp += chop;
        BT += chop;
    }

    //  Compute basic fields.
    target->sign = a->sign == b->sign; // Sign is positive if signs are equal.
    target->exp  = Aexp + Bexp;        // Add the exponents.

    //  Allocate mantissa
    if(target->coef == a->coef) {
      // "In-place" mul
      if(AL+BL > target->len) {
        target->coef = (uint32_t*) realloc(target->coef, sizeof(uint32_t)*(AL+BL));
        // Prevent use-after-free
        if(AT == BT) {
          BT = target->coef;
        }
        AT = target->coef;
      }
    } else if(target->coef == b->coef) {
      if(AL+BL > target->len) {
        target->coef = (uint32_t*) realloc(target->coef, sizeof(uint32_t)*(AL+BL));
        // Prevent use-after-free
        if(AT == BT) {
          AT = target->coef;
        }
        BT = target->coef;
      }
    } else {
      if(target->coef != NULL) {
        free(target->coef);
      }
      target->coef = (uint32_t*) malloc(sizeof(uint32_t)*(AL+BL));
    }

    target->len  = AL + BL; // Add the lenghts for now. May need to correct later.

    uint32_t *CT = target->coef;
    uint32_t CL = target->len;

    if(CL == 2) {
      // Fast path for really small multiplications
      uint64_t r = (uint64_t)AT[0]*BT[0];
      CT[0] = r & 0xffffffff;
      CT[1] = r >> 32;
    } else if(CL < BASECASE_THRESH) {
      _basecase_mul(CT, CL, AT, AL, BT, BL);
    } else {
      uint32_t AC = checksum(AT, AL);
      uint32_t BC = checksum(BT, BL);
      uint32_t CC_should = checksum_mul(AC, BC);
      int bits = _fft_mul(CT, CL, AT, AL, BT, BL);
      if(CC_should != checksum(CT, CL)) {
        fprintf(stderr, "Multiplication error (CL=%d in %d bits)\n", CL, bits);
        abort();
      }
    }
    #ifdef DEBUG
    bigfloat_print("t", target);
    #endif

    //  Check top word and correct length.
    if (target->coef[target->len - 1] == 0) {
      target->len--;
    }
}
