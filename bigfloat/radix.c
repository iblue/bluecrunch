#include <malloc.h>
#include <string.h> // memcpy
#include <stdlib.h> // abort
#include "math.h"
#include "bigfloat.h"

#define OLDBASE 4294967296
#define NEWBASE 100000000
#define LOG (log(OLDBASE)/log(NEWBASE))

// Calculates approximate number of digits in base 100000000 of a number.
size_t bigfloat_radix_decimals(bigfloat_t target) {
  double approx = 0.0;

  if(target->len >= 3) {
    approx += (double)target->coef[target->len-3];
  }

  if(target->len >= 2) {
    approx += (double)target->coef[target->len-2] * (double)OLDBASE;
  }

  approx += (double)target->coef[target->len-1] * (double)OLDBASE * (double)OLDBASE;

  double digits = log(approx)/log(NEWBASE);
  double additional = ((double)target->len-3.0)*log(OLDBASE)/log(NEWBASE);

  double epsilon = 1e-10;

  return ceil(digits + additional + epsilon);
}

// XXX untested
int64_t _scale(bigfloat_t a) {
  // Already scaled
  if(a->exp == 0) {
    return 0;
  }

  double scale = -LOG*a->exp;
  int64_t intscale = scale;

  //printf("bigfloat_radix: %f =~ %ld\n", scale, intscale);

  // Multiply by scaling.
  bigfloat_t mul;
  bigfloat_new(mul);
  bigfloat_set(mul, 100000000);
  bigfloat_exp(mul, mul, intscale, 8);
  bigfloat_mul(a, a, mul, 0, 8);

  // FIXME: Check if exp is zero
  bigfloat_floor(a);

  return -intscale;
}

void _radix_recurse(uint32_t *coef, size_t len) {
  if(len <= 1) {
    return;
  }

  bigfloat_t high;
  bigfloat_new(high);

  bigfloat_t low;
  bigfloat_new(low);

  // We are cheating here...
  bigfloat_t src;
  src->coef = coef;
  src->len  = len;
  src->exp  = 0;
  src->sign = 1;

  size_t split = 1;
  while(split < len) {
    split <<= 1;
  }
  split >>= 1;

  //printf("%ld split to (%ld,%ld)\n", len, split, len-split);

  // TODO: Performance: This will be in the precalculated base conv table
  bigfloat_t exp;
  bigfloat_new(exp);
  bigfloat_set(exp, NEWBASE);
  bigfloat_exp(exp, exp, split, 8);

  // FIXME: Performance: We calculate more than we need.
  // FIXME: Performance: Precalculate bigfloat_rcp(exp)
  bigfloat_div(high, src, exp, 0, 8);
  bigfloat_floor(high);

  bigfloat_mul(low, exp, high, 0, 8);
  bigfloat_neg(low);
  bigfloat_add(low, low, src, 0);

  if(low->len + high->len != len) {
    fprintf(stderr, "radix conversion error\n");
    abort();
  }

  bigfloat_realloc(low, bigfloat_radix_decimals(low));
  bigfloat_realloc(high, bigfloat_radix_decimals(high));

  // Write target
  memcpy(coef, low->coef, low->len*sizeof(uint32_t));
  memcpy(coef+low->len, high->coef, high->len*sizeof(uint32_t));

  size_t lowlen = low->len;
  size_t highlen = high->len;

  bigfloat_free(low);
  bigfloat_free(high);

  // Recurse
  _radix_recurse(coef, lowlen);
  _radix_recurse(coef+lowlen, highlen);
}

// Converts from OLDBASE to NEWBASE
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  if(target->coef != a->coef) {
    bigfloat_copy(target, a);
  }

  // Scale digits
  int64_t newexp = _scale(target);

  size_t new_digits = bigfloat_radix_decimals(target);
  bigfloat_realloc(target, new_digits);

  // This is the recursion
  _radix_recurse(target->coef, new_digits);

  if(newexp != 0) {
    target->exp = newexp;
  }
}
