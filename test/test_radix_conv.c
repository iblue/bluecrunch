#include <malloc.h>
#include <assert.h>
#include "bigfloat.h"

int main(void) {
  {
    // = 111222333444555666777888999000
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = 0;
    a->sign = 1;
    a->len  = 4;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0xab1f4658;
    a->coef[1] = 0x79875a70;
    a->coef[2] = 0x6760f539;
    a->coef[3] = 0x1;

    //bigfloat_print("a", a);

    bigfloat_t b;
    bigfloat_new(b);

    bigfloat_radix(b, a);
    //bigfloat_print10("b", b);

    assert(b->len  == 4);
    assert(b->sign == 1);
    assert(b->exp  == 0);
    assert(b->coef[0] == 88999000);
    assert(b->coef[1] == 56667778);
    assert(b->coef[2] == 33344455);
    assert(b->coef[3] == 111222);
  }

  return 0;
}
