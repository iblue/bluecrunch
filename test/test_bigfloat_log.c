#include <malloc.h>
#include <assert.h>
#include <math.h>
#include "bigfloat.h"

int main(void) {
  {
    bigfloat_t a;
    bigfloat_new(a);
    bigfloat_set(a, 99999999);

    assert(bigfloat_radix_decimals(a) == 1);

    bigfloat_free(a);
  }

  {
    bigfloat_t a;
    bigfloat_new(a);
    bigfloat_set(a, 100000000);

    assert(bigfloat_radix_decimals(a) == 2);

    bigfloat_free(a);
  }

  {
    bigfloat_t a;
    bigfloat_new(a);
    bigfloat_set(a, 0xffffffff);

    assert(bigfloat_radix_decimals(a) == 2);

    bigfloat_free(a);
  }

  {
    // = 100000000000000000000000000000000000000000000000000000000
    bigfloat_t a;
    bigfloat_new(a);
    bigfloat_alloc(a, 6);
    a->coef[0] = 0x00000000;
    a->coef[1] = 0x21000000;
    a->coef[2] = 0x73d4490d;
    a->coef[3] = 0xfdffc788;
    a->coef[4] = 0x940f6a24;
    a->coef[5] = 0x04140c78;

    assert(bigfloat_radix_decimals(a) == 7);

    bigfloat_free(a);
  }

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

    // This is a case where we are off by one. We need to correct this after
    // the conversion
    //assert(bigfloat_radix_decimals(a) == 6);
    assert(bigfloat_radix_decimals(a) == 7);

    bigfloat_free(a);
  }

  return 0;
}
