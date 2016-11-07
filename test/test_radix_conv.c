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

  {
    // 607851afc94c9f559.e395724f752a9800
    // = 111222333444555666777.888999123
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = -2;
    a->sign = 1;
    a->len  = 5;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0x752aa200;
    a->coef[1] = 0xe395724f;
    a->coef[2] = 0x94c9f559;
    a->coef[3] = 0x07851afc;
    a->coef[4] = 0x6;

    //bigfloat_print("a", a);

    bigfloat_t b;
    bigfloat_new(b);

    bigfloat_radix(b, a);
    //bigfloat_print10("b", b);

    assert(b->len  == 5);
    assert(b->sign == 1);
    assert(b->exp  == -2);
    assert(b->coef[0] == 30000000);
    assert(b->coef[1] == 88899912);
    assert(b->coef[2] == 55666777);
    assert(b->coef[3] == 23334445);
    assert(b->coef[4] == 11122);
  }

  {
    // This is e (to ~100 digits)
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = -13;
    a->sign = 1;
    a->len  = 14;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0]  = 0xda06c806;
    a->coef[1]  = 0x8c31d763;
    a->coef[2]  = 0xf4bf8d8d;
    a->coef[3]  = 0x926cfbe5;
    a->coef[4]  = 0x324e7738;
    a->coef[5]  = 0x5190cfef;
    a->coef[6]  = 0xa784d904;
    a->coef[7]  = 0x38b4da56;
    a->coef[8]  = 0x62e7160f;
    a->coef[9]  = 0x9cf4f3c7;
    a->coef[10] = 0xbf715880;
    a->coef[11] = 0x8aed2a6a;
    a->coef[12] = 0xb7e15162;
    a->coef[13] = 0x00000002;

    bigfloat_radix(a, a);
    bigfloat_print10("a", a);
  }


  return 0;
}
