#include <malloc.h>
#include "bigfloat.h"

#define CEIL_POS(X) ((X-(int)(X)) > 0 ? (int)(X+1) : (int)(X))
#define CEIL_NEG(X) ((X-(int)(X)) < 0 ? (int)(X-1) : (int)(X))
#define CEIL(X) ( ((X) > 0) ? CEIL_POS(X) : CEIL_NEG(X) )

// Converts to base 10
// (Ok, we are converting from base 4294967296 to 100000000)
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  // This is log_100000000(4294967296)
  //double scale = 1.204119982655924780854955578897972107072759525848434165241;
  // Now multiply by 1000000000^exponent
  //assert(exponent == 0); // Fix that later.

  // FIXME: Need bigfloat_exp first.
  size_t limbs = 2; // FIXME!
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
}
