#include <malloc.h>
#include <string.h> // memmove
#include "math.h"
#include "bigfloat.h"

#define OLDBASE 4294967296
#define NEWBASE 100000000
#define LOG (log(OLDBASE)/log(NEWBASE))

// Calculates approximate number of digits in base 100000000 of a number.
size_t bigfloat_radix_decimals(bigfloat_t target) {
  double approx = 0.0;

  if(target->len >= 3) {
    approx += (double)target->coef[target->len-3] * pow((double)OLDBASE, target->len-3);
  }

  if(target->len >= 2) {
    approx += (double)target->coef[target->len-2] * pow((double)OLDBASE, target->len-2);
  }

  approx += (double)target->coef[target->len-1] * pow((double)OLDBASE, target->len-1);

  double digits = log(approx)/log(NEWBASE);

  double epsilon = 1e-10;

  return ceil(digits + epsilon);
}

// Converts from OLDBASE to NEWBASE
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  bigfloat_t scaled;
  bigfloat_new(scaled);
  int64_t newexp = 0;


  if(newexp != 0) {
    target->exp = newexp;

    // Fix overflow by carry.
    if(target->coef[target->len-1] > NEWBASE) {
      target->len++;
      target->coef = (uint32_t*) realloc(target->coef, target->len*sizeof(uint32_t));
      target->coef[target->len-1] = target->coef[target->len-2]/NEWBASE;
      target->coef[target->len-2] %= NEWBASE;
    }
  }
}
