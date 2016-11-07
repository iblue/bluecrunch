#include <malloc.h>
#include "bigfloat.h"

void bigfloat_alloc(bigfloat_t target, size_t size) {
  bigfloat_free(target);

  target->exp     = 0;
  target->sign    = 1;
  target->len     = size;
  target->coef    = (uint32_t*) malloc(target->len*sizeof(uint32_t));
}
