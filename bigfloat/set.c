#include <malloc.h>
#include "bigfloat.h"

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
