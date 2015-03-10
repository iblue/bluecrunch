//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
#include "bigfloat.h"
#include "fft.h"

void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds) {
    //  Multiplication

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

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
    target->sign = a->sign == b->sign;  //  Sign is positive if signs are equal.
    target->exp  = Aexp + Bexp;       //  Add the exponents.
    target->len    = AL + BL;           //  Add the lenghts for now. May need to correct later.

    //  Allocate mantissa
    if(target->coef != NULL) {
      free(target->coef);
    }
    target->coef = (uint32_t*)malloc(sizeof(uint32_t)*(target->len));

    if(target->len == 2) {
      #ifdef DDEBUG
      printf("fastpath\n");
      #endif
      // Fast path for really small multiplications
      uint64_t result = (uint64_t)AT[0]*BT[0];
      target->coef[0] = result % 1000000000;
      result /= 1000000000;
      target->coef[1] = result;
    } else if(target->len < 350) {
      #ifdef DDEBUG
      printf("basecase\n");
      #endif
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

          carry  = value / 1000000000;
          value %= 1000000000;
          target->coef[i+j] = value;
        }
        for(size_t j=i+BL;j<target->len;j++) {
          if(carry == 0) {
            break;
          }
          value = target->coef[j];
          value += carry;

          carry  = value / 1000000000;
          value %= 1000000000;
          target->coef[j] = value;
        }
      }
    } else {
      #ifdef DDEBUG
      printf("fft\n");
      #endif
      //  Perform multiplication.
      int digits_per_point;
      if(target->len > 80000000) {
        digits_per_point = 2;
      } else if(target->len > 1000000) {
        digits_per_point = 3;
      } else {
        digits_per_point = 4;
      }

      int points_per_word = 9/digits_per_point;
      if(9%digits_per_point) {
        points_per_word++;
      }

      //  Determine minimum FFT size.
      int k = 0;
      size_t length = 1;
      while (length < points_per_word*target->len) {
        length <<= 1;
        k++;
      }

      //  Allocate FFT arrays
      __m128d *Ta = (__m128d*)_mm_malloc(length * sizeof(__m128d), 32);
      __m128d *Tb = (__m128d*)_mm_malloc(length * sizeof(__m128d), 32);

      //  Make sure the twiddle table is big enough.

      int_to_fft(Ta,k,AT,AL, digits_per_point);           //  Convert 1st operand
      int_to_fft(Tb,k,BT,BL, digits_per_point);           //  Convert 2nd operand
      if (twiddle_table_size - 1 < k) {
        fft_forward_uncached(Ta,k,tds);
        fft_forward_uncached(Tb,k,tds);
      } else {
        fft_forward(Ta,k,tds);
        fft_forward(Tb,k,tds);
      }
      fft_pointwise(Ta,Tb,k);//  Pointwise multiply
      _mm_free(Tb);
      if (twiddle_table_size - 1 < k) {
        fft_inverse_uncached(Ta,k,tds);
      } else {
        fft_inverse(Ta,k,tds);
      }
      fft_to_int(Ta,k,target->coef,target->len, digits_per_point);   //  Convert back to word array.
      _mm_free(Ta);
    }
    #ifdef DDEBUG
    for(size_t i=AL;i-->0;) {
      printf("%09d", AT[i]);
    }
    printf(",");
    for(size_t i=BL;i-->0;) {
      printf("%09d", BT[i]);
    }
    printf(",");
    for(size_t i=target->len;i-->0;) {
      printf("%09d", target->coef[i]);
    }
    printf("\n");
    #endif

    //  Check top word and correct length.
    if (target->coef[target->len - 1] == 0) {
      target->len--;
    }
}
