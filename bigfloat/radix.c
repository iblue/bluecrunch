#include <malloc.h>
#include "bigfloat.h"

// Converts to base 10
// (Ok, we are converting from base 4294967296 to 1000000000)
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  // This is -log_1000000000(4294967296)
  double scale = -1.07032887347193313853773829235375298406467513408749703577;
  uint64_t exponent = scale*a->exp; // floor.

  // Now multiply by 1000000000^exponent
  //assert(exponent == 0); // Fix that later.

  // FIXME: Need bigfloat_exp first.
}

// Converts.
size_t inline bigfloat_radix_recursion(uint32_t *target, size_t len) {
  
}
