#include <malloc.h>
#include <assert.h>
#include "bigfloat.h"

int main(void) {
  /*
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
  */

  {
    // This is e.
    bigfloat_t a;
    bigfloat_new(a);
    a->exp  = -22;
    a->sign = 1;
    a->len  = 23;
    a->coef = (uint32_t*) malloc(a->len*sizeof(uint32_t));
    a->coef[0] = 0xfd5f24d6;
    a->coef[1] = 0xd55c4d79;
    a->coef[2] = 0xf7b46bce;
    a->coef[3] = 0x158d9554;
    a->coef[4] = 0x7c19bb42;
    a->coef[5] = 0x90cfd47d;
    a->coef[6] = 0x57f59584;
    a->coef[7] = 0x4f7c7b57;
    a->coef[8] = 0xbb1185eb;
    a->coef[9] = 0xda06c80a;
    a->coef[10] = 0x8c31d763;
    a->coef[11] = 0xf4bf8d8d;
    a->coef[12] = 0x926cfbe5;
    a->coef[13] = 0x324e7738;
    a->coef[14] = 0x5190cfef;
    a->coef[15] = 0xa784d904;
    a->coef[16] = 0x38b4da56;
    a->coef[17] = 0x62e7160f;
    a->coef[18] = 0x9cf4f3c7;
    a->coef[19] = 0xbf715880;
    a->coef[20] = 0x8aed2a6a;
    a->coef[21] = 0xb7e15162;
    a->coef[22] = 0x00000002;

    bigfloat_t b;
    bigfloat_new(b);

    bigfloat_radix(b, a, 8);

    assert(b->len  == 27);
    assert(b->sign == 1);
    assert(b->exp  == -26);

    // FIXME: Fix length
    assert(b->coef[0] == 15738341);
    assert(b->coef[1] == 25101901);
    assert(b->coef[2] == 80753195);
    assert(b->coef[3] == 32338298);
    assert(b->coef[4] == 94349076);
    assert(b->coef[5] == 32328627);
    assert(b->coef[6] == 56307381);
    assert(b->coef[7] == 29526059);
    assert(b->coef[8] == 72900334);
    assert(b->coef[9] == 66290435);
    assert(b->coef[10] == 81741359);
    assert(b->coef[11] == 3059921);
    assert(b->coef[12] == 63919320);
    assert(b->coef[13] == 42742746);
    assert(b->coef[14] == 78525166);
    assert(b->coef[15] == 45713821);
    assert(b->coef[16] == 35354759);
    assert(b->coef[17] == 24076630);
    assert(b->coef[18] == 69676277);
    assert(b->coef[19] == 95957496);
    assert(b->coef[20] == 47093699);
    assert(b->coef[21] == 24977572);
    assert(b->coef[22] == 47135266);
    assert(b->coef[23] == 35360287);
    assert(b->coef[24] == 84590452);
    assert(b->coef[25] == 71828182);
    assert(b->coef[26] == 2);
  }

  /*
  {
    // = 99999999999999999999999999999999999999999999999999999999
    bigfloat_t a;
    bigfloat_new(a);
    bigfloat_alloc(a, 6);
    a->coef[0] = 0xffffffff;
    a->coef[1] = 0x20ffffff;
    a->coef[2] = 0x73d4490d;
    a->coef[3] = 0xfdffc788;
    a->coef[4] = 0x940f6a24;
    a->coef[5] = 0x04140c78;

    bigfloat_radix(a, a);

    assert(a->len  == 7);
    assert(a->sign == 1);
    assert(a->exp  == 0);
    assert(a->coef[0] == 99999999);
    assert(a->coef[1] == 99999999);
    assert(a->coef[2] == 99999999);
    assert(a->coef[3] == 99999999);
    assert(a->coef[4] == 99999999);
    assert(a->coef[5] == 99999999);
    assert(a->coef[6] == 99999999);

    bigfloat_free(a);
  }
  */


  return 0;
}
