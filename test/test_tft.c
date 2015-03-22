#include "fft.h"
#include <assert.h>
#include <stdio.h>

static inline void assert_fp(complex double a, complex double b) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi)\n",
      creal(b), cimag(b), creal(a), cimag(b));
    abort();
  }
}


int main() {
  fft_ensure_table(8); // up to 256 values

  {
    __attribute__ ((aligned (32))) complex double values[] = {6, 2, -2+2*I, 0};

    tft_inverse(values, 3);

    assert_fp(values[0], 4);
    assert_fp(values[1], 8);
    assert_fp(values[2], 12);
    assert_fp(values[3], 0);
  }

  /*
   * Forward FFT:
   *
   * stage 0   1.000+0.000i   2.000+0.000i   3.000+0.000i   4.000+0.000i   5.000+0.000i  6.000+0.000i  7.000+0.000i   0.000+0.000i
   * stage 1   6.000+0.000i   8.000+0.000i  10.000+0.000i   4.000+0.000i  -4.000+0.000i -2.828-2.828i -0.000-4.000i  -2.828+2.828i
   * stage 2  16.000+0.000i  12.000+0.000i  -4.000+0.000i   0.000+4.000i  -4.000-4.000i -5.657+0.000i -4.000+4.000i   5.657-0.000i
   * stage 3  28.000+0.000i   4.000+0.000i  -4.000+4.000i  -4.000-4.000i  -9.657-4.000i  1.657-4.000i  1.657+4.000i  -9.657+4.000i
   */
  { 
    __attribute__ ((aligned (32))) complex double values[] = {
      28,
      4,
      -4 + 4*I,
      -4 - 4*I,
      -9.6568542494923797 - 4*I,
       1.6568542494923797 - 4*I,
       1.6568542494923801 + 4*I,
      0, //-9.6568542494923797 + 4*I
    };


    tft_inverse(values, 7);

    assert_fp(values[0], 8);
    assert_fp(values[1], 16);
    assert_fp(values[2], 24);
    assert_fp(values[3], 32);
    assert_fp(values[4], 40);
    assert_fp(values[5], 48);
    assert_fp(values[6], 56);
    assert_fp(values[7], 0);
  }

  /*
   * Forward FFT:
   *
   * stage 0:  2.000+0i  11.000+0i  17.000+00i   7.000+00i   5.000+00.000i 13.000+00.000i  0.000+00.000i   0.000+00.000i
   * stage 1:  7.000+0i  24.000+0i  17.000+00i   7.000+00i  -3.000+00.000i -1.414-01.414i  0.000+17.000i  -4.950+04.950i
   * stage 2: 24.000+0i  31.000+0i -10.000+00i   0.000+17i  -3.000+17.000i -6.364+03.536i -3.000-17.000i   6.364+03.536i
   * stage 3: 55.000+0i  -7.000+0i -10.000+17i -10.000-17i  -9.364+20.536i  3.364+13.464i  3.364-13.464i  -9.364-20.536i
   */
  {
    __attribute__ ((aligned (32))) complex double values[] = {
       55,
      -7
      -10+17*I,
      -10-17*I,
      -9.36396103+20.53553391*I,
       3.36396103+13.46446609*I,
      0, // 3.36396103-13.46446609*I,
      0, //-9.36396103-20.53553391*I,
    };

    tft_inverse(values, 6);

    assert_fp(values[0],  7*8);
    assert_fp(values[1], 11*8);
    assert_fp(values[2], 17*8);
    assert_fp(values[3],  7*8);
    assert_fp(values[4],  5*8);
    assert_fp(values[5], 13*8);
    assert_fp(values[6], 0);
    assert_fp(values[7], 0);
  }
}
