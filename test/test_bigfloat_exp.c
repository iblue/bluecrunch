#include <malloc.h>
#include <assert.h>
#include "bigfloat.h"

int main(void) {
  {
    //2^256
    bigfloat_t a, b;
    bigfloat_new(a);
    bigfloat_set(a, 0x1234);

    bigfloat_new(b);

    bigfloat_exp(b, a, 0xd7);

    // Expected result: 0x10298d253fecb89d491933590a085405544b42f1d9c94
    // 0ca666a0b296a5c661fbcac60e9bbb478e26b3e640d61c2832cd6b61f7d67ea0
    // 2df7c14ce8279e958b612d26d54ae6e5a7ce6e7c4b24156ff511648e672df968
    // 070750ec664af4ed6b56a11935c75d52781051282e03e2b55c4e9f9eef659e32
    // 7a3ad9e7c1711345d7e26c9db89d3eb94e59f53f2c8924737736cbd1321f8d80
    // ea626aaa0a8b9691847816ebfa322982ebd48bca0475d049403e81e57362d491
    // b6625b25de3b364fe1ee0c4ddf80121e042a81231e258d627343b35c08b2e51e
    // aa632b034342209cb69c7171de9c9bdc47e90d39de0b73703f46d15c8cc1bf19
    // 14405dc7de708e4c8814342828b0ba1f04c4cc854a672040e4a93c9400000000
    // 0000000000000000000000000000000000000000000000000000000000000000
    // 00000000000000000000000000000000000
    assert(b->len  == 82);
    assert(b->sign == 1);
    assert(b->exp  == 0);
    assert(b->coef[81] == 0x10298d25);
    assert(b->coef[80] == 0x3fecb89d);
    assert(b->coef[79] == 0x49193359);
    assert(b->coef[78] == 0x0a085405);
    assert(b->coef[77] == 0x544b42f1);
    assert(b->coef[76] == 0xd9c940ca);
    assert(b->coef[75] == 0x666a0b29);
    assert(b->coef[74] == 0x6a5c661f);
    assert(b->coef[73] == 0xbcac60e9);
    assert(b->coef[72] == 0xbbb478e2);
    assert(b->coef[71] == 0x6b3e640d);
    assert(b->coef[70] == 0x61c2832c);
    assert(b->coef[69] == 0xd6b61f7d);
    assert(b->coef[68] == 0x67ea02df);
    assert(b->coef[67] == 0x7c14ce82);
    assert(b->coef[66] == 0x79e958b6);
    assert(b->coef[65] == 0x12d26d54);
    assert(b->coef[64] == 0xae6e5a7c);
    assert(b->coef[63] == 0xe6e7c4b2);
    assert(b->coef[62] == 0x4156ff51);
    assert(b->coef[61] == 0x1648e672);
    assert(b->coef[60] == 0xdf968070);
    assert(b->coef[59] == 0x750ec664);
    assert(b->coef[58] == 0xaf4ed6b5);
    assert(b->coef[57] == 0x6a11935c);
    assert(b->coef[56] == 0x75d52781);
    assert(b->coef[55] == 0x051282e0);
    assert(b->coef[54] == 0x3e2b55c4);
    assert(b->coef[53] == 0xe9f9eef6);
    assert(b->coef[52] == 0x59e327a3);
    assert(b->coef[51] == 0xad9e7c17);
    assert(b->coef[50] == 0x11345d7e);
    assert(b->coef[49] == 0x26c9db89);
    assert(b->coef[48] == 0xd3eb94e5);
    assert(b->coef[47] == 0x9f53f2c8);
    assert(b->coef[46] == 0x92473773);
    assert(b->coef[45] == 0x6cbd1321);
    assert(b->coef[44] == 0xf8d80ea6);
    assert(b->coef[43] == 0x26aaa0a8);
    assert(b->coef[42] == 0xb9691847);
    assert(b->coef[41] == 0x816ebfa3);
    assert(b->coef[40] == 0x22982ebd);
    assert(b->coef[39] == 0x48bca047);
    assert(b->coef[38] == 0x5d049403);
    assert(b->coef[37] == 0xe81e5736);
    assert(b->coef[36] == 0x2d491b66);
    assert(b->coef[35] == 0x25b25de3);
    assert(b->coef[34] == 0xb364fe1e);
    assert(b->coef[33] == 0xe0c4ddf8);
    assert(b->coef[32] == 0x0121e042);
    assert(b->coef[31] == 0xa81231e2);
    assert(b->coef[30] == 0x58d62734);
    assert(b->coef[29] == 0x3b35c08b);
    assert(b->coef[28] == 0x2e51eaa6);
    assert(b->coef[27] == 0x32b03434);
    assert(b->coef[26] == 0x2209cb69);
    assert(b->coef[25] == 0xc7171de9);
    assert(b->coef[24] == 0xc9bdc47e);
    assert(b->coef[23] == 0x90d39de0);
    assert(b->coef[22] == 0xb73703f4);
    assert(b->coef[21] == 0x6d15c8cc);
    assert(b->coef[20] == 0x1bf19144);
    assert(b->coef[19] == 0x05dc7de7);
    assert(b->coef[18] == 0x08e4c881);
    assert(b->coef[17] == 0x4342828b);
    assert(b->coef[16] == 0x0ba1f04c);
    assert(b->coef[15] == 0x4cc854a6);
    assert(b->coef[14] == 0x72040e4a);
    assert(b->coef[13] == 0x93c94000);
    assert(b->coef[12] == 0x00000000);
    assert(b->coef[11] == 0x00000000);
    assert(b->coef[10] == 0x00000000);
    assert(b->coef[ 9] == 0x00000000);
    assert(b->coef[ 8] == 0x00000000);
    assert(b->coef[ 7] == 0x00000000);
    assert(b->coef[ 6] == 0x00000000);
    assert(b->coef[ 5] == 0x00000000);
    assert(b->coef[ 4] == 0x00000000);
    assert(b->coef[ 3] == 0x00000000);
    assert(b->coef[ 2] == 0x00000000);
    assert(b->coef[ 1] == 0x00000000);
    assert(b->coef[ 0] == 0x00000000);

    bigfloat_free(a);
    bigfloat_free(b);
  }

  return 0;
}
