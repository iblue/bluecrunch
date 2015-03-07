#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <iostream>

//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
extern "C" {
  #include "fft.h"
}
#include "bigfloat.h"

void bigfloat_new(BigFloat target) {
  target->sign = true;
  target->exp  = 0;
  target->L    = 0;
  target->T    = NULL;
}

void bigfloat_set(BigFloat target, uint32_t x, bool sign_) {
  bigfloat_free(target);

  target->exp  = 0;

  if (x == 0) {
      target->L    = 0;
      target->sign = true;
      target->T    = NULL;
      return;
  }

  target->sign = sign_;

  target->T = (uint32_t*) malloc(sizeof(uint32_t));
  target->T[0] = x;
  target->L    = 1;
}

void bigfloat_free(BigFloat target) {
  if(target->T != NULL) {
    free(target->T);
    target->T = NULL;
  }
}

////////////////////////////////////////////////////////////////////////////////
//  String Conversion
char* bigfloat_to_string_trimmed(const BigFloat value, size_t digits, int64_t* ret_exp) {
    //  Converts this object to a string with "digits" significant figures.

    //  After calling this function, the following expression is equal to the
    //  numeric value of this object. (after truncation of precision)
    //      str + " * 10^" + (return value)

    char* str;

    // Handle zero
    if (value->L == 0) {
      *ret_exp = 0;
      str = (char*) malloc(2);
      str[0] = '0';
      str[1] = '\0';
      return str;
    }

    //  Collect operands
    int64_t exponent = value->exp;
    size_t length = value->L;
    uint32_t *ptr = value->T;

    if (digits == 0){
        //  Use all digits.
        digits = length * 9;
    } else {
        //  Truncate precision
        size_t words = (digits + 17) / 9;
        if (words < length){
            size_t chop = length - words;
            exponent += chop;
            length = words;
            ptr += chop;
        }
    }
    exponent *= 9;

    //  Build string
    size_t c = length;
    str = (char*)malloc(9*length+1);
    while (c-- > 0){
      uint32_t word = ptr[c];

      for (int i = 8; i >= 0; i--){
          str[9*c+i] = word % 10 + '0';
          word /= 10;
      }
    }

    //  Count leading zeros
    size_t leading_zeros = 0;
    while (str[leading_zeros] == '0') {
      leading_zeros++;
    }
    digits += leading_zeros;

    //  Truncate
    if (digits < strlen(str)){
      exponent += strlen(str) - digits;
      str[digits] = '\0';
    }

    *ret_exp = exponent;
    return str;
}

size_t bigfloat_to_string(char* string, const BigFloat value, size_t digits) {
    //  Convert this number to a string. Auto-select format type.
    if (value->L == 0) {
      string = (char*)malloc(2);
      string[0] = '0';
      string[1] = '\0';
      return 2;
    }

    int64_t mag = value->exp + value->L;

    //  Use scientific notation of out of range.
    if (mag > 1 || mag < 0) {
      fprintf(stderr, "Scientific notation required\n");
      abort();
    }

    //  Convert
    int64_t exponent;
    char * str = bigfloat_to_string_trimmed(value, digits, &exponent);

    //  Less than 1
    /*
    if (mag == 0){
      std::string ret;
      if (value->sign) {
        ret = std::string("0.") + str;
      } else {
        ret = std::string("-0.") + str;
      }

      string = (char*) ret.c_str();
      return ret.size();
    }
    */
    if (mag == 0) {
      fprintf(stderr, "Not implemented\n");
      abort();
    }

    //  Get a string with the digits before the decimal place.
    char before_decimal[10];
    sprintf(before_decimal, "%d", value->T[value->L - 1]);

    //  Nothing after the decimal place.
    if (exponent >= 0){
      /*
      std::string ret;
        if (value->sign){
            ret = before_decimal + ".";
        }else{
            ret = std::string("-") + before_decimal + ".";
        }

      string = (char*) ret.c_str();
      return ret.size();
      */
      fprintf(stderr, "Not implemented\n");
      abort();
    }

    //  Get digits after the decimal place.
    std::string after_decimal = str.substr((size_t)(strlen(str) + exponent),(size_t)-exponent);

    std::string ret;
    if (value->sign){
        ret = before_decimal + "." + after_decimal;
    }else{
        ret = std::string("-") + before_decimal + "." + after_decimal;
    }

    string = (char*) ret.c_str();
    return ret.size();
}

////////////////////////////////////////////////////////////////////////////////
//  Getters
uint32_t _bigfloat_word_at(const BigFloat target, int64_t mag) {
  //  Returns the word at the mag'th digit place.
  //  This is useful for additions where you need to access a specific "digit place"
  //  of the operand without having to worry if it's out-of-bounds.

  //  This function is mathematically equal to:
  //      (return value) = floor(this * (10^9)^-mag) % 10^9

  if (mag < target->exp)
      return 0;
  if (mag >= target->exp + (int64_t)target->L)
      return 0;
  return target->T[(size_t)(mag - target->exp)];
}

int _bigfloat_ucmp(const BigFloat a, const BigFloat b) {
    //  Compare function that ignores the sign.
    //  This is needed to determine which direction subtractions will go.

    //  Magnitude
    int64_t magA = a->exp + a->L;
    int64_t magB = b->exp + b->L;
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
void bigfloat_negate(BigFloat num) {
  if(num->L == 0) {
    return;
  }

  num->sign = !num->sign;
}

void _bigfloat_uadd(BigFloat target, const BigFloat a, const BigFloat b, size_t p) {
    //  Perform addition ignoring the sign of the two operands.

    //  Magnitude
    int64_t magA = a->exp + a->L;
    int64_t magB = b->exp + b->L;
    int64_t top = std::max(magA, magB);
    int64_t bot = std::min(a->exp, b->exp);

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
    target->L     = (uint32_t)TL;

    //  Allocate mantissa
    target->T = (uint32_t*) malloc(sizeof(uint32_t)*(TL + 1));

    //  Add
    uint32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++){
      uint32_t word = _bigfloat_word_at(a, bot) + _bigfloat_word_at(b, bot) + carry;
      carry = 0;
      if (word >= 1000000000){
          word -= 1000000000;
          carry = 1;
      }
      target->T[c] = word;
    }

    //  Carry out
    if (carry != 0) {
        target->T[target->L++] = 1;
    }
}

void _bigfloat_usub(BigFloat target, const BigFloat a, const BigFloat b, size_t p) {
    //  Perform subtraction ignoring the sign of the two operands.

    //  "this" must be greater than or equal to x. Otherwise, the behavior
    //  is undefined.

    //  Magnitude
    int64_t magA = a->exp + a->L;
    int64_t magB = b->exp + b->L;
    int64_t top = std::max(magA, magB);
    int64_t bot = std::min(a->exp, b->exp);

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
    target->L     = (uint32_t)TL;

    //  Allocate mantissa
    target->T = (uint32_t*) malloc(sizeof(uint32_t)*(TL + 1));

    //  Subtract
    int32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++) {
        int32_t word = (int32_t)_bigfloat_word_at(a, bot) - (int32_t)_bigfloat_word_at(b, bot) - carry;
        carry = 0;
        if (word < 0){
            word += 1000000000;
            carry = 1;
        }
        target->T[c] = word;
    }

    //  Strip leading zeros
    while (target->L > 0 && target->T[target->L - 1] == 0)
        target->L--;
    if (target->L == 0){
        target->exp = 0;
        target->sign = true;
        free(target->T);
        target->T = NULL;
    }
}

void bigfloat_add(BigFloat target, const BigFloat a, const BigFloat b, size_t p) {
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
        bigfloat_negate(target);
      }
    }
}

void bigfloat_sub(BigFloat target, const BigFloat a, const BigFloat b, size_t p) {
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
      bigfloat_negate(target);
    }
  }
}

void bigfloat_mul(BigFloat target, const BigFloat a, const BigFloat b, size_t p, int tds) {
    //  Multiplication

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Either operand is zero.
    if (a->L == 0 || b->L == 0) {
      bigfloat_new(target);
      return;
    }

    if (p == 0) {
        //  Default value. No trunction.
        p = a->L + b->L;
    } else {
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    //  Collect operands.
    int64_t Aexp = a->exp;
    int64_t Bexp = b->exp;
    size_t AL = a->L;
    size_t BL = b->L;
    uint32_t *AT = a->T;
    uint32_t *BT = b->T;

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
    target->L    = AL + BL;           //  Add the lenghts for now. May need to correct later.

    //  Allocate mantissa
    if(target->T != NULL) {
      free(target->T);
    }
    target->T = (uint32_t*)malloc(sizeof(uint32_t)*(target->L));

    //  Perform multiplication.

    //  Determine minimum FFT size.
    int k = 0;
    size_t length = 1;
    while (length < 5*target->L) {
        length <<= 1;
        k++;
    }

    //  Perform a convolution using FFT.
    //  Yeah, this is slow for small sizes, but it's asympotically optimal.

    // FIXME: Comment incorrect. 2 digits per point!
    //  3 digits per point is small enough to not encounter round-off error
    //  until a transform size of 2^30.
    //  A transform length of 2^29 allows for the maximum product size to be
    //  2^29 * 3 = 1,610,612,736 decimal digits.
    /*if (k > 29)
        throw "FFT size limit exceeded.";*/

    //  Allocate FFT arrays
    __m128d *Ta = (__m128d*)_mm_malloc(length * sizeof(__m128d), 16);
    __m128d *Tb = (__m128d*)_mm_malloc(length * sizeof(__m128d), 16);

    //  Make sure the twiddle table is big enough.
    if (twiddle_table_size - 1 < k) {
      throw "Table is not large enough.";
    }

    int_to_fft(Ta,k,AT,AL);           //  Convert 1st operand
    int_to_fft(Tb,k,BT,BL);           //  Convert 2nd operand
    fft_forward(Ta,k,tds); //  Transform 1st operand
    fft_forward(Tb,k,tds); //  Transform 2nd operand
    fft_pointwise(Ta,Tb,k);//  Pointwise multiply
    fft_inverse(Ta,k,tds); //  Perform inverse transform.
    fft_to_int(Ta,k,target->T,target->L);   //  Convert back to word array.
    _mm_free(Ta);
    _mm_free(Tb);

    //  Check top word and correct length.
    if (target->T[target->L - 1] == 0)
        target->L--;
}

void bigfloat_rcp(BigFloat target, const BigFloat a, size_t p, int tds) {
    //  Compute reciprocal using Newton's Method.

    //  r1 = r0 - (r0 * x - 1) * r0

    if (a->L == 0) {
      fprintf(stderr, "Divide by Zero\n");
      abort();
    }

    //  Collect operand
    int64_t Aexp = a->exp;
    size_t AL = a->L;
    uint32_t *AT = a->T;

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

        //  Rebuild a BigFloat.
        uint64_t val64 = (uint64_t)val;

        target->sign = a->sign;

        target->T = (uint32_t*)malloc(sizeof(uint32_t)*2);
        target->T[0] = (uint32_t)(val64 % 1000000000);
        target->T[1] = (uint32_t)(val64 / 1000000000);
        target->L = 2;
        target->exp = Aexp;

        return;
    }

    //  Half the precision
    size_t s = p / 2 + 1;
    if (p == 1) s = 0;
    if (p == 2) s = 1;

    //  Recurse at half the precision
    BigFloat T;
    bigfloat_new(T);
    bigfloat_rcp(T, a, s, tds);

    //  r1 = r0 - (r0 * x - 1) * r0
    BigFloat one;
    bigfloat_new(one);
    bigfloat_set(one, 1, 1);
    BigFloat tmp;
    bigfloat_new(tmp);
    bigfloat_mul(tmp, a, T, p, tds);
    BigFloat tmp2;
    bigfloat_new(tmp2);
    bigfloat_sub(tmp2, tmp, one, p);
    BigFloat tmp3;
    bigfloat_new(tmp3);
    bigfloat_mul(tmp3, tmp2, T, p, tds);
    bigfloat_sub(target, T, tmp3, p);

    bigfloat_free(tmp);
    bigfloat_free(tmp2);
    bigfloat_free(tmp3);
    bigfloat_free(T);
    bigfloat_free(one);
}

void bigfloat_div(BigFloat target, const BigFloat a, const BigFloat b, size_t p, int tds) {
  //  Division
  BigFloat rcp;
  bigfloat_new(rcp);
  bigfloat_rcp(rcp, b, p, tds);
  bigfloat_mul(target, a, rcp, p, tds);
  bigfloat_free(rcp);
}
