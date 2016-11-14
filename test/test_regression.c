#include <assert.h>
#include <malloc.h>
#include "bigfloat.h"
#include "fft.h"


int main() {
  fft_ensure_table(32);

  /*
  {
    bigfloat_t a;

    a->exp = 0;
    a->len = 4;
    a->sign = 1;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0x44210000;
    a->coef[1] = 0xb1514c28;
    a->coef[2] = 0xce88deb8;
    a->coef[3] = 0x00000002;

    bigfloat_t target;
    bigfloat_new(target);

    bigfloat_mul(target, a, a, 0, 8);

    assert(target->exp == 0);
    assert(target->len == 7);
    assert(target->sign == 1);
    assert(target->coef[0] == 0x00000000);
    assert(target->coef[1] == 0xf4718c41);
    assert(target->coef[2] == 0x21c8adaa);
    assert(target->coef[3] == 0x28b6543c);
    assert(target->coef[4] == 0x68a5465b);
    assert(target->coef[5] == 0xe0c40a81);
    assert(target->coef[6] == 0x7);
  }
  */
  {
    bigfloat_t a;

    a->exp = 0;
    a->len = 4;
    a->sign = 1;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0x44210000;
    a->coef[1] = 0xb1514c28;
    a->coef[2] = 0xce88deb8;
    a->coef[3] = 0x00000002;


    bigfloat_mul(a, a, a, 0);

    assert(a->exp == 0);
    assert(a->len == 7);
    assert(a->sign == 1);
    assert(a->coef[0] == 0x00000000);
    assert(a->coef[1] == 0xf4718c41);
    assert(a->coef[2] == 0x21c8adaa);
    assert(a->coef[3] == 0x28b6543c);
    assert(a->coef[4] == 0x68a5465b);
    assert(a->coef[5] == 0xe0c40a81);
    assert(a->coef[6] == 0x7);

    bigfloat_free(a);
  }
}
