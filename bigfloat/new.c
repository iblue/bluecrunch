#include <malloc.h>
#include "bigfloat.h"

void bigfloat_new(bigfloat_t target) {
  target->sign = 1;
  target->exp  = 0;
  target->len  = 0;
  target->coef = NULL;
}
