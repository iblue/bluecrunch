#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>

//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
#include "bigfloat.h"

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

// Gets word at a given magnitude if set.
uint32_t _bigfloat_word_at(const bigfloat_t target, int64_t mag) {
  if (mag < target->exp) {
    return 0;
  }

  if (mag >= target->exp + (int64_t)target->len) {
    return 0;
  }

  return target->coef[(size_t)(mag - target->exp)];
}

// Compares the absolute values of a and b. Returns:
//  1 if the absolute of a is bigger
// -1 if the absolute of b is bigger
//  0 if the absolutes are equal
int _bigfloat_ucmp(const bigfloat_t a, const bigfloat_t b) {
  int64_t a_mag = a->exp + a->len;
  int64_t b_mag = b->exp + b->len;

  if (a_mag > b_mag) {
    return 1;
  }

  if (a_mag < b_mag) {
    return -1;
  }

  int64_t mag = a_mag;

  while(mag >= a->exp || mag >= b->exp) {
    uint32_t a_word = _bigfloat_word_at(a, mag);
    uint32_t b_word = _bigfloat_word_at(b, mag);

    if (a_word < b_word) {
      return -1;
    }

    if (a_word > b_word) {
      return 1;
    }

    mag--;
  }

  return 0;
}

void _bigfloat_uadd(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p) {
  int64_t a_mag = a->exp + a->len;
  int64_t b_mag = b->exp + b->len;
  int64_t top   = max(a_mag, b_mag);
  int64_t bot   = min(a->exp, b->exp);

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
    int64_t a_mag = a->exp + a->len;
    int64_t b_mag = b->exp + b->len;
    int64_t top = max(a_mag, b_mag);
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

