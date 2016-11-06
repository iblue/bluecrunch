#include <assert.h>
#include "bigfloat.h"


int main() {
  {
    bigfloat_t a;
    bigfloat_t target;

    a->exp = 0;
    a->len = 4;
    a->sign = 1;
    a->coef[0] = 0x44210000;
    a->coef[1] = 0xb1514c28;
    a->coef[2] = 0xce88deb8;
    a->coef[3] = 0x00000002;


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
}
