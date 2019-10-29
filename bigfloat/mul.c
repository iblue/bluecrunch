//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>
#include <stdlib.h>

#include <cilk/cilk.h>
#include "bigfloat.h"
#include "fft.h"

#include <math.h>
#include <alloca.h>

#define BASECASE_THRESH 100
#define KARATZUBA_THRESH 3300

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

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

uint32_t checksum_add(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p  = (av+bv)%CHECKPRIME;
  return (uint32_t)p;
}

uint32_t checksum_sub(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p;
  if(av > bv) {
    p = av-bv;
  } else {
    p = av+CHECKPRIME-bv;
  }
  p  = p%CHECKPRIME;
  return (uint32_t)p;
}

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
int _fft_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
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

void _basecase_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
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

  if(AL == 0 || BL == 0 || CL == 0) {
    return;
  }

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

void debug_print(const char* name, uint32_t* T, size_t L) {
  printf("%s = COEF {len = %ld, coef = ", name, L);
  printf("{ ");
  for(size_t i=0; i<L; i++) {
    printf("0x%08x", T[i]);
    if(i != L - 1) {
      printf(", ");
    }
  }
  printf(" }}\n");
}

void static inline _ll_add(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
#ifdef DEBUG
  uint32_t CA = checksum(AT, AL);
  uint32_t CB = checksum(BT, BL);
#endif

  uint32_t carry = 0;
  for (size_t i = 0; i<CL; i++) {
    uint64_t word = (uint64_t)carry;
    if(i<AL) {
      word += (uint64_t)AT[i];
    }
    if(i<BL) {
      word += (uint64_t)BT[i];
    }
    carry = 0;
    if (word >= 0x100000000){
      word -= 0x100000000;
      carry = 1;
    }
    CT[i] = word;
  }

#ifdef DEBUG
  if(checksum(CT, CL) != checksum_add(CA, CB)) {
    printf("Checksum Error in _ll_add\n");
    debug_print("A", AT, AL);
    debug_print("B", BT, BL);
    debug_print("C", CT, CL);
    abort();
  }
#endif
}

void static inline _ll_sub_inplace(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL) {
#ifdef DEBUG
  uint32_t CCb = checksum(CT, CL);
  uint32_t CA  = checksum(AT, AL);
#endif

  int32_t carry = 0;
  for (size_t i = 0; i<CL; i++) {
    int64_t word = (int64_t)CT[i]-(int64_t)carry;
    if(i<AL) {
      word -= (int64_t)AT[i];
    }
    carry = 0;
    if (word < 0){
      word += 0x100000000;
      carry = 1;
    }
    CT[i] = word;
  }

#ifdef DEBUG
  if(checksum(CT, CL) != checksum_sub(CCb, CA)) {
    printf("Checksum Error in _ll_sub_inplace\n");
    debug_print("C", CT, CL);
    abort();
  }
#endif
}

void _karatzuba_mul(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL, uint32_t *BT, size_t BL) {
  if(CL < BASECASE_THRESH || AL == 0 || BL == 0 || CL == 0) {
    _basecase_mul(CT, CL, AT, AL, BT, BL);
    return;
  }

#ifdef DEBUG
  uint32_t AC = checksum(AT, AL);
  uint32_t BC = checksum(BT, BL);
  uint32_t CC_should = checksum_mul(AC, BC);
#endif

  char inplaceOverride = (CT == AT || CT == BT);

  if(inplaceOverride) {
    uint32_t* PT = malloc(CL*sizeof(uint32_t));
    _karatzuba_mul(PT, CL, AT, AL, BT, BL);
    memcpy(CT, PT, CL*sizeof(uint32_t));
    free(PT);
    return;
  }

  // Zero target (FIXME: Needed?)
  //memset(CT, 0, sizeof(uint32_t)*CL);

  // Karatzuba algorithm:
  //
  // 1. Split A and B in two halfes: (AHi, ALo), (BHi, BLo).
  size_t SL = max(AL, BL)/2; // This is the split length.

  uint32_t* AHi  = AT+SL;
  size_t    AHiL = AL>SL ? AL-SL : 0;
  uint32_t* ALo  = AT;
  size_t    ALoL = AL>SL ? SL : AL;

  uint32_t* BHi  = BT+SL;
  size_t    BHiL = BL>SL ? BL-SL : 0;
  uint32_t* BLo  = BT;
  size_t    BLoL = BL>SL ? SL : BL;

#ifdef DEBUG
  printf("DEBUG: Karatzuba AL=%ld, BL=%ld, CL=%ld, SL=%ld\n", AL, BL, CL, SL);

  debug_print("A",   AT,  AL);
  debug_print("AHi", AHi, AHiL);
  debug_print("ALo", ALo, ALoL);

  debug_print("B",   BT,  BL);
  debug_print("BHi", BHi, BHiL);
  debug_print("BLo", BLo, BLoL);
#endif

  // Target also composed of two halfes, we call them Z0, Z2;
  uint32_t* Z2  = CT+SL*2;
  size_t    Z2L = CL-SL*2; //assert(CL>=SL*2);
  uint32_t* Z0  = CT;
  size_t    Z0L = SL*2;


  // Compute:
  // Z2 = AHi*BHi
  // Z0 = ALo*BLo
  // Z1 = (AHi+ALo)*(BHi+BLo) - Z2 - Z0
  // Result = Z2 * B^(2*SL) + Z1 * B^(SL) + Z0
  //
  // We save Z2 and Z0 directly in CT with some index shifts
  // Then we compute Z1 in scrach space and add it with some shifts.

  // Z2 = AHi*BHi. Save in upper half of CL
  _karatzuba_mul(Z2, Z2L, AHi, AHiL, BHi, BHiL);
#ifdef DEBUG
  debug_print("Z2 (AHi*BHi)", Z2, Z2L);
#endif
  // Z0 = ALo*BLo. Save in lower half of CL
  _karatzuba_mul(Z0, Z0L, ALo, ALoL, BLo, BLoL);
#ifdef DEBUG
  debug_print("Z0 (ALo*BLo)", Z0, Z0L);
#endif

// We allocate the buffers in once piece in production mode, because malloc is expensive.
#ifdef DEBUG
  size_t    T0L = max(AHiL, ALoL)+1;  // AHi+ALo (AHi: len AL-SL, ALo: len SL, SL<AL-SL)
  uint32_t* T0 = malloc(T0L*sizeof(uint32_t));
  size_t    T1L = max(BHiL, BLoL)+1; // equiv to A
  uint32_t* T1 = malloc(T1L*sizeof(uint32_t));
  size_t    T2L = T0L+T1L;  // will contain T0*T1
  uint32_t* T2 = malloc(T2L*sizeof(uint32_t));
#else
  size_t    T0L = max(AHiL, ALoL)+1;  // AHi+ALo (AHi: len AL-SL, ALo: len SL, SL<AL-SL)
  size_t    T1L = max(BHiL, BLoL)+1; // equiv to A
  size_t    T2L = T0L+T1L;  // will contain T0*T1
  uint32_t* T = alloca((T0L+T1L+T2L)*sizeof(uint32_t));
  uint32_t* T0 = T;
  uint32_t* T1 = T0+T0L;
  uint32_t* T2 = T1+T1L;
#endif

#ifdef DEBUG
  printf("DEBUG: lens T0=%ld, T1=%ld, T2=%ld\n", T0L, T1L, T2L);
#endif

  // T0 = AHi+ALo
  _ll_add(T0, T0L, AHi, AHiL, ALo, ALoL);

  // T1 = BHi+BLo
  _ll_add(T1, T1L, BHi, BHiL, BLo, BLoL);

  // T2 = T0*T1
  _karatzuba_mul(T2, T2L, T0, T0L, T1, T1L);

#ifdef DEBUG
  debug_print("T2 (T0*T1)", T2, T2L);
#endif

  // T2 -= Z0
#ifdef DEBUG
  printf("T2L = %ld, Z0L = %ld\n", T2L, SL*2);
#endif
  _ll_sub_inplace(T2, T2L, Z0, Z0L);

  // T2 -= Z2
  // FIXME: Optimize. If AL<SL || BL<SL, AHi = 0 || BHi = 0, therefore Z2 = 0
  // and this can be skipped. Leave in for debugging for now.
  _ll_sub_inplace(T2, T2L, Z2, Z2L);

  // CT += T2*Base^SL (Base^SL is implemented via shift of CT index)
  {
    uint32_t carry = 0;
    for (size_t i = 0; i<CL-SL; i++){
      uint64_t word = (uint64_t)CT[i+SL] + (uint64_t)carry;
      if(i<T2L) {
        word += (uint64_t)T2[i];
      }
      carry = 0;
      if (word >= 0x100000000){
        word -= 0x100000000;
        carry = 1;
      }
      CT[i+SL] = word;
    }
  }

#ifdef DEBUG
  debug_print("C", CT, CL);
#endif

#ifdef DEBUG
  free(T2);
  free(T1);
  free(T0);
#else
  //alloca allocates on stack
  //free(T);
#endif

#ifdef DEBUG
  if(CC_should != checksum(CT, CL)) {
    fprintf(stderr, "Karatzuba Internal Multiplication error (CL=%ld)\n", CL);
    abort();
  }
#endif
}

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
    } else if(CL < KARATZUBA_THRESH) {
      _karatzuba_mul(CT, CL, AT, AL, BT, BL);
    } else {
      uint32_t AC = checksum(AT, AL);
      uint32_t BC = checksum(BT, BL);
      uint32_t CC_should = checksum_mul(AC, BC);
      int bits = _fft_mul(CT, CL, AT, AL, BT, BL);
      if(CC_should != checksum(CT, CL)) {
        fprintf(stderr, "FFT Multiplication error (CL=%d in %d bits)\n", CL, bits);
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
