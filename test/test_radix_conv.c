#include <malloc.h>
#include "bigfloat.h"

int main(void) {
  {
    // = 123456789123456789
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = 0;
    a->sign = 1;
    a->len  = 2;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0xacd05f15;
    a->coef[1] = 0x01b69b4b;

    bigfloat_print("a", a);

    bigfloat_t b;
    bigfloat_new(b);

    bigfloat_radix(b, a);
    bigfloat_print("b", b);
    // expected: 0x499602d2 (twice)
  }

  return 0;
}
