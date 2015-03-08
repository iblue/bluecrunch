#define _USE_MATH_DEFINES
#include <malloc.h>
#include "bigfloat.h"

void bigfloat_zero(bigfloat_t target) {
  target->exp   = 0;
  target->len   = 0;
  target->sign  = 1;
  if(target->coef) {
    free(target->coef);
  }
  target->coef  = NULL;
}

int bigfloat_iszero(bigfloat_t target) {
  return target->len == 0;
}
