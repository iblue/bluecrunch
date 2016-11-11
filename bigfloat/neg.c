#include "bigfloat.h"

void bigfloat_neg(bigfloat_t target) {
  if(bigfloat_iszero(target)) {
    return;
  }

  target->sign = !target->sign;
}
