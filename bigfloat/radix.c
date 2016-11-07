#include <malloc.h>
#include <string.h> // memmove
#include "math.h"
#include "bigfloat.h"

// Converts to base 10
// (Ok, we are converting from base 4294967296 to 100000000)
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  // This is log_100000000(4294967296)
  double log = 1.204119982655924780854955578897972107072759525848434165241;

  bigfloat_t scaled;
  bigfloat_new(scaled);
  int64_t newexp = 0;

  if(a->exp != 0) {
    double scale = -log*a->exp;
    uint64_t intscale = scale; // Off by one errors!
    newexp = -scale;

    // Multiply by scaling.
    bigfloat_t mul;
    bigfloat_new(mul);
    bigfloat_set(mul, 100000000);
    bigfloat_exp(mul, mul, intscale, 8);
    bigfloat_mul(scaled, mul, a, 0, 8);

    // FIXME: Check if exp is zero
    bigfloat_floor(scaled);
  } else {
    // FIXME: Just link
    bigfloat_copy(scaled, a);
  }

  // Actually that could be refined to (may require 1 limb less)
  // size_t limbs = ceil(log(a->coef[a->len-1])/log + (a->len-1)*log)
  size_t limbs = ceil(scaled->len*log);
  size_t upper_len = limbs/2;

  // Convert upper
  bigfloat_t exp;
  bigfloat_new(exp);
  bigfloat_set(exp, 100000000);
  bigfloat_exp(exp, exp, upper_len, 8);

  bigfloat_t high;
  bigfloat_new(high);
  // FIXME: Performance: We calculate more than we need for the floor.
  bigfloat_div(high, scaled, exp, upper_len, 8); // Contains now upper half
  bigfloat_floor(high);

  // Convert lower
  bigfloat_t low;
  bigfloat_new(low);
  bigfloat_mul(low, exp, high, 0, 8);
  bigfloat_neg(low);
  bigfloat_add(low, low, scaled, 0);

  // Recurse
  if(low->len > 1) {
    bigfloat_radix(low, low);
  }

  if(high->len > 1) {
    bigfloat_radix(high, high);
  }

  // Write target
  bigfloat_alloc(target, low->len + high->len);
  memmove(target->coef, low->coef, low->len*sizeof(uint32_t));
  memmove(target->coef+low->len, high->coef, high->len*sizeof(uint32_t));

  // FIXME: We can skip this by writing directly into target by manipulating
  // some shit.
  bigfloat_free(low);
  bigfloat_free(high);

  target->exp = newexp;
}
