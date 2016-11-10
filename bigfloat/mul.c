//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
#include "bigfloat.h"
#include "fft.h"

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

// Multiplies C = A*B using FFT where C = CL integers in CT, etc.
// CT, CL = target ptr and length
// AT, AL = source 1 ptr and length
// BT, BL = source 2 ptr and length
// tds = numbers of threads to be used
void static inline _fft_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL, int tds) {
  int bits_per_point = 12;
  int points_per_word = 3;

  //  Determine minimum FFT size.
  int k = 0;
  size_t length = 1;
  while (length < points_per_word*CL) {
    length <<= 1;
    k++;
  }

  if(length/2 > points_per_word*CL+10) {
    length /= 2;
  }

  // If the arguments are the same, we skip one conversion.
  char needB = (AT != BT || AL != BL);

  //  Allocate FFT arrays
  complex double *Ta = (complex double*)_mm_malloc(length * sizeof(complex double), 32);
  complex double *Tb = NULL;

  if(needB) {
    Tb = (complex double*)_mm_malloc(length * sizeof(complex double), 32);
  }

  //  Convert Numbers to FFT
  int_to_fft(Ta, k, AT, AL, bits_per_point);

  if(needB) {
    int_to_fft(Tb, k, BT, BL, bits_per_point);
  }

  // If numbers are too big, use multiplication without table.
  if (twiddle_table_size - 1 < k) {
    fft_forward_uncached(Ta, k, tds);
    if(needB) {
      fft_forward_uncached(Tb, k, tds);
    }
  } else {
    fft_forward(Ta, k, tds);
    if(needB) {
      fft_forward(Tb, k, tds);
    }
  }

  // Pointwise multiply
  if(needB) {
    fft_pointwise(Ta, Tb, k);
  } else {
    fft_pointwise(Ta, Ta, k);
  }

  // Inverse transform
  if (twiddle_table_size - 1 < k) {
    fft_inverse_uncached(Ta, k, tds);
  } else {
    fft_inverse(Ta, k, tds);
  }

  // Convert including carryout
  fft_to_int(Ta, k, CT, CL, bits_per_point);

  // Free FFT arrays
  if(needB) {
    _mm_free(Tb);
  }
  _mm_free(Ta);
}

void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds) {
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

    /*#ifndef DEBUG
    #define BASECASE_OPTIMIZATION
    #endif*/
    #ifdef BASECASE_OPTIMIZATION
    if(target->len == 2) {
      // Fast path for really small multiplications
      uint64_t result = (uint64_t)AT[0]*BT[0];
      target->coef[0] = result & 0xffffffff;
      target->coef[1] = result >> 32;
    } else if(target->len < 350 && target->coef != a->coef && target->coef != b->coef) {
      // FIXME: This breaks if we are multiplying in-place. Fix that.
      for(size_t i=0;i<target->len;i++) {
        target->coef[i] = 0;
      }
      for(size_t i=0;i<AL;i++) {
        uint64_t carry = 0;
        uint64_t value;
        for(size_t j=0;j<BL;j++) {
          value = target->coef[i+j];
          value += (uint64_t) AT[i] * (uint64_t) BT[j];
          value += carry;

          carry  = value >> 32;
          value &= 0xffffffff;
          target->coef[i+j] = value;
        }
        for(size_t j=i+BL;j<target->len;j++) {
          if(carry == 0) {
            break;
          }
          value = target->coef[j];
          value += carry;

          carry  = value / UINT32_MAX;
          value %= UINT32_MAX;
          target->coef[j] = value;
        }
      }
    } else
    #endif
    {
      uint32_t *CT = target->coef;
      uint32_t CL = target->len;

      _fft_mul(CT, CL, AT, AL, BT, BL, tds);
    }
    #ifdef DEBUG
    bigfloat_print("t", target);
    #endif

    //  Check top word and correct length.
    if (target->coef[target->len - 1] == 0) {
      target->len--;
    }
}
