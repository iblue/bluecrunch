#include <math.h>
#include <stdint.h>
#include <string.h>
#include <malloc.h>
#include "bigfloat.h"


void bigfloat_floor(bigfloat_t target) {
  if(target->coef == NULL) {
    return;
  }

  if(target->exp >= 0) {
    return;
  }

  size_t shift = -target->exp;
  target->len -= shift;
  target->exp = 0;
  memmove(target->coef, target->coef+shift, sizeof(uint32_t)*target->len);
  target->coef = (uint32_t*) realloc(target->coef, sizeof(uint32_t)*target->len);
}
