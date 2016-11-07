#include <malloc.h>
#include <string.h> // memmove
#include "math.h"
#include "bigfloat.h"

// Converts to base 10
// (Ok, we are converting from base 4294967296 to 100000000)
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  // This is log_100000000(4294967296)
  double scale = 1.204119982655924780854955578897972107072759525848434165241;
  // Now multiply by 1000000000^exponent
  //assert(exponent == 0); // Fix that later.

  // Actually that could be refined to (may require 1 limb less)
  // size_t limbs = ceil(log(a->coef[a->len-1])/scale + (a->len-1)*scale)
  size_t limbs = ceil(a->len*scale);
  size_t upper_len = limbs/2;

  // Convert upper
  bigfloat_t exp;
  bigfloat_new(exp);
  bigfloat_set(exp, 100000000);
  bigfloat_exp(exp, exp, upper_len, 8);

  bigfloat_t high;
  bigfloat_new(high);
  // FIXME: Performance: We calculate more than we need for the floor.
  bigfloat_div(high, a, exp, upper_len, 8); // Contains now upper half
  bigfloat_floor(high);

  // Convert lower
  bigfloat_t low;
  bigfloat_new(low);
  bigfloat_mul(low, exp, high, 0, 8);
  bigfloat_neg(low);
  bigfloat_add(low, low, a, 0);

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
}
