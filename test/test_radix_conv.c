#include <malloc.h>
#include "bigfloat.h"

int main(void) {
  {
    //2^256
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = 0;
    a->sign = 1;
    a->len  = 8;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 1;
    a->coef[1] = 0;
    a->coef[2] = 0;
    a->coef[3] = 0;
    a->coef[4] = 0;
    a->coef[5] = 0;
    a->coef[6] = 0;
    a->coef[7] = 0;

    bigfloat_print("a", a);
  }

  return 0;
}
