#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>

//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
#include "fft.h"
#include "bigfloat.h"

uint32_t _bigfloat_word_at(const bigfloat_t target, int64_t mag);
int _bigfloat_ucmp(const bigfloat_t a, const bigfloat_t b);

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

void bigfloat_new(bigfloat_t target) {
  target->sign = 1;
  target->exp  = 0;
  target->len  = 0;
  target->coef  = NULL;
}

void bigfloat_neg(bigfloat_t num) {
  if(bigfloat_iszero(num)) {
    return;
  }

  num->sign = !num->sign;
}

void bigfloat_set(bigfloat_t target, uint32_t value) {
  bigfloat_free(target);

  if (value == 0) {
    bigfloat_zero(target);
    return;
  }

  target->exp     = 0;
  target->sign    = 1;
  target->len     = 1;
  target->coef    = (uint32_t*) malloc(target->len*sizeof(uint32_t));
  target->coef[0] = value;
}

void bigfloat_free(bigfloat_t target) {
  if(target->coef != NULL) {
    free(target->coef);
    target->coef = NULL;
  }
}

size_t int_to_str(uint32_t val, char* str) {
  for (int i=8;i>=0;i--){
    str[i] = val % 10 + '0';
    val   /= 10;
  }
  return 9;
}

// Skipps leading zeros
size_t int_to_str_trimmed(uint32_t val, char* str) {
  if(val == 2) {
    str[0] = '2';
    return 1;
  } else {
    fprintf(stderr, "Not implemented\n");
    abort();
  }
}

// Returns length of string, fills char* string with value
size_t bigfloat_to_string(char* string, const bigfloat_t value, size_t digits) {
  char* initial_string = string;
  if(value->len == 0) {
    string[0] = '0';
    string[1] = '\0';
    return 1;
  }

  int mag = value->len + value->exp;

  size_t c = value->len-1;

  if(mag == 1) {
    string += int_to_str_trimmed(value->coef[c], string);
  }

  *string++ = '.';

  size_t min_c;

  if(value->len > digits/9) {
    min_c = value->len-digits/9-1;
  } else {
    min_c = 0;
  }

  while(c-->min_c) {
    string += int_to_str(value->coef[c], string);
  }

  *string++ = '\0';

  // Truncate if required
  if(string - initial_string > digits) {
    initial_string[digits+1] = '\0';
    return digits;
  }

  return string - initial_string-1;
}

////////////////////////////////////////////////////////////////////////////////
//  Getters
uint32_t _bigfloat_word_at(const bigfloat_t target, int64_t mag) {
  //  Returns the word at the mag'th digit place.
  //  This is useful for additions where you need to access a specific "digit place"
  //  of the operand without having to worry if it's out-of-bounds.

  //  This function is mathematically equal to:
  //      (return value) = floor(this * (10^9)^-mag) % 10^9

  if (mag < target->exp)
      return 0;
  if (mag >= target->exp + (int64_t)target->len)
      return 0;
  return target->coef[(size_t)(mag - target->exp)];
}

int _bigfloat_ucmp(const bigfloat_t a, const bigfloat_t b) {
    //  Compare function that ignores the sign.
    //  This is needed to determine which direction subtractions will go.

    //  Magnitude
    int64_t magA = a->exp + a->len;
    int64_t magB = b->exp + b->len;
    if (magA > magB)
        return 1;
    if (magA < magB)
        return -1;

    //  Compare
    int64_t mag = magA;
    while (mag >= a->exp || mag >= b->exp){
        uint32_t wordA = _bigfloat_word_at(a, mag);
        uint32_t wordB = _bigfloat_word_at(b, mag);
        if (wordA < wordB)
            return -1;
        if (wordA > wordB)
            return 1;
        mag--;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  Arithmetic
void _bigfloat_uadd(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
    //  Perform addition ignoring the sign of the two operands.

    //  Magnitude
    int64_t magA = a->exp + a->len;
    int64_t magB = b->exp + b->len;
    int64_t top = max(magA, magB);
    int64_t bot = min(a->exp, b->exp);

    //  Target length
    int64_t TL = top - bot;

    if (p == 0) {
        //  Default value. No trunction.
        p = (size_t)TL;
    } else {
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    //  Perform precision truncation.
    if (TL > (int64_t)p){
        bot = top - p;
        TL = p;
    }

    //  Compute basic fields.
    target->sign  = a->sign;
    target->exp   = bot;
    target->len     = (uint32_t)TL;

    //  Allocate mantissa
    target->coef = (uint32_t*) malloc(sizeof(uint32_t)*(TL + 1));

    //  Add
    uint32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++){
      uint32_t word = _bigfloat_word_at(a, bot) + _bigfloat_word_at(b, bot) + carry;
      carry = 0;
      if (word >= 1000000000){
          word -= 1000000000;
          carry = 1;
      }
      target->coef[c] = word;
    }

    //  Carry out
    if (carry != 0) {
        target->coef[target->len++] = 1;
    }
}

void _bigfloat_usub(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
    //  Perform subtraction ignoring the sign of the two operands.

    //  "this" must be greater than or equal to x. Otherwise, the behavior
    //  is undefined.

    //  Magnitude
    int64_t magA = a->exp + a->len;
    int64_t magB = b->exp + b->len;
    int64_t top = max(magA, magB);
    int64_t bot = min(a->exp, b->exp);

    //  Truncate precision
    int64_t TL = top - bot;

    if (p == 0) {
        //  Default value. No trunction.
        p = (size_t)TL;
    } else {
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    if (TL > (int64_t)p){
        bot = top - p;
        TL = p;
    }

    //  Compute basic fields.
    target->sign  = a->sign;
    target->exp   = bot;
    target->len     = (uint32_t)TL;

    //  Allocate mantissa
    target->coef = (uint32_t*) malloc(sizeof(uint32_t)*(TL + 1));

    //  Subtract
    int32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++) {
        int32_t word = (int32_t)_bigfloat_word_at(a, bot) - (int32_t)_bigfloat_word_at(b, bot) - carry;
        carry = 0;
        if (word < 0){
            word += 1000000000;
            carry = 1;
        }
        target->coef[c] = word;
    }

    //  Strip leading zeros
    while (target->len > 0 && target->coef[target->len - 1] == 0) {
      target->len--;
    }

    if (target->len == 0){
      bigfloat_zero(target);
    }
}

void bigfloat_add(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
    //  Addition

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Same sign. Add.
    if (a->sign == b->sign) {
      _bigfloat_uadd(target, a, b, p);
    } else { // Differing signs. Subtract.
      //  this > x
      if (_bigfloat_ucmp(a, b) > 0) {
        _bigfloat_usub(target, a, b, p);
      } else { //  this < x
        _bigfloat_usub(target, b, a, p);
        bigfloat_neg(target);
      }
    }
}

void bigfloat_sub(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
  //  Subtraction

  //  The target precision is p.
  //  If (p = 0), then no truncation is done. The entire operation is done
  //  at maximum precision with no data loss.

  //  Different sign. Add.
  if (a->sign != b->sign) {
    _bigfloat_uadd(target, a, b, p);
  } else { // Differing signs. Subtract.
    //  this > x
    if (_bigfloat_ucmp(a, b) > 0) {
      _bigfloat_usub(target, a, b, p);
    } else { //  this < x
      _bigfloat_usub(target, b, a, p);
      bigfloat_neg(target);
    }
  }
}

void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds) {
    //  Multiplication

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Either operand is zero.
    if (a->len == 0 || b->len == 0) {
      bigfloat_new(target);
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
      #ifdef DEBUG
      printf("fastpath\n");
      #endif
      // Fast path for really small multiplications
      uint64_t result = (uint64_t)AT[0]*BT[0];
      target->coef[0] = result % 1000000000;
      result /= 1000000000;
      target->coef[1] = result;
    } else if(target->len < 500) {
      #ifdef DEBUG
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
      #ifdef DEBUG
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
      __m128d *Ta = (__m128d*)_mm_malloc(length * sizeof(__m128d), 16);
      __m128d *Tb = (__m128d*)_mm_malloc(length * sizeof(__m128d), 16);

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
    #ifdef DEBUG
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

void bigfloat_rcp(bigfloat_t target, const bigfloat_t a, size_t p, int tds) {
    //  Compute reciprocal using Newton's Method.

    //  r1 = r0 - (r0 * x - 1) * r0

    if (a->len == 0) {
      fprintf(stderr, "Divide by Zero\n");
      abort();
    }

    //  Collect operand
    int64_t Aexp = a->exp;
    size_t AL = a->len;
    uint32_t *AT = a->coef;

    //  End of recursion. Generate starting point.
    if (p == 0){
        //  Truncate precision to 3.
        p = 3;
        if (AL > p){
            size_t chop = AL - p;
            AL = p;
            Aexp += chop;
            AT += chop;
        }

        //  Convert number to floating-point.
        double val = AT[0];
        if (AL >= 2)
            val += AT[1] * 1000000000.;
        if (AL >= 3)
            val += AT[2] * 1000000000000000000.;

        //  Compute reciprocal.
        val = 1. / val;
        Aexp = -Aexp;

        //  Scale
        while (val < 1000000000.){
            val *= 1000000000.;
            Aexp--;
        }

        //  Rebuild a bigfloat_t.
        uint64_t val64 = (uint64_t)val;

        target->sign = a->sign;

        target->coef = (uint32_t*)malloc(sizeof(uint32_t)*2);
        target->coef[0] = (uint32_t)(val64 % 1000000000);
        target->coef[1] = (uint32_t)(val64 / 1000000000);
        target->len = 2;
        target->exp = Aexp;

        return;
    }

    //  Half the precision
    size_t s = p / 2 + 1;
    if (p == 1) s = 0;
    if (p == 2) s = 1;

    //  Recurse at half the precision
    bigfloat_t T;
    bigfloat_new(T);
    bigfloat_rcp(T, a, s, tds);

    //  r1 = r0 - (r0 * x - 1) * r0
    bigfloat_t one;
    bigfloat_new(one);
    bigfloat_set(one, 1);
    bigfloat_t tmp;
    bigfloat_new(tmp);
    bigfloat_mul(tmp, a, T, p, tds);
    bigfloat_t tmp2;
    bigfloat_new(tmp2);
    bigfloat_sub(tmp2, tmp, one, p);
    bigfloat_free(tmp);
    bigfloat_t tmp3;
    bigfloat_new(tmp3);
    bigfloat_mul(tmp3, tmp2, T, p, tds);
    bigfloat_free(tmp2);
    bigfloat_sub(target, T, tmp3, p);
    bigfloat_free(tmp3);
    bigfloat_free(T);
    bigfloat_free(one);
}

void bigfloat_div(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds) {
  //  Division
  bigfloat_t rcp;
  bigfloat_new(rcp);
  bigfloat_rcp(rcp, b, p, tds);
  bigfloat_mul(target, a, rcp, p, tds);
  bigfloat_free(rcp);
}