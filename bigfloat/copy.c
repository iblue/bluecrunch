#include <malloc.h>
#include <string.h>
#include "bigfloat.h"

void bigfloat_copy(bigfloat_t target, const bigfloat_t source) {
  bigfloat_free(target);

  target->exp  = source->exp;
  target->sign = source->sign;
  target->len  = source->len;
  target->coef = (uint32_t*) malloc(target->len*sizeof(uint32_t));
  memcpy(target->coef, source->coef, sizeof(uint32_t)*target->len);
}
