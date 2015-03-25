#include "fft.h"
#include <assert.h>
#include <stdio.h>

#define assert_fp(a, b) _assert_fp(a, b, __FILE__, __LINE__)

static inline void _assert_fp(complex double a, complex double b, char* file, int line) {
  double radius = 1e-6;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d\n",
      creal(b), cimag(b), creal(a), cimag(b), file, line);
    abort();
  }
}


int main() {
  fft_ensure_table(8); // up to 256 values

  {
    __attribute__ ((aligned (32))) complex double values[] = {1, 2, 3, 0};

    tft_forward(values, 3, 2);

    assert_fp(values[0], 6);
    assert_fp(values[1], 2);
    assert_fp(values[2], -2+2*I);
    assert_fp(values[3], 0);
  }


  {
    __attribute__ ((aligned (32))) complex double values[] = {6, 2, -2+2*I, 0};

    tft_inverse(values, 3, 2);

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


    tft_inverse(values, 7, 3);

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
      -7,
      -10+17*I,
      -10-17*I,
      -9.36396103+20.53553391*I,
       3.36396103+13.46446609*I,
      0, // 3.36396103-13.46446609*I,
      0, //-9.36396103-20.53553391*I,
    };

    tft_inverse(values, 6, 3);

    assert_fp(values[0],  2*8);
    assert_fp(values[1], 11*8);
    assert_fp(values[2], 17*8);
    assert_fp(values[3],  7*8);
    assert_fp(values[4],  5*8);
    assert_fp(values[5], 13*8);
    assert_fp(values[6], 0);
    assert_fp(values[7], 0);
  }

  /*
   * Forward FFT:
   *
   * stage 0:    1  34   26      95      |  53                  0.000+0.000i     0.000+00.000i     00.000+00.000i
   * stage 1:   54  34   26      95      | -52                 24.042+24.042i    0.000+26.000i    -67.175+67.175i
   * stage 2:   80  129  28         -61i | -52    +26i        -43.134+91.217i  -52.000-26.000i     43.134+91.217i
   * stage 3:  209 -49   28-61i  28+61i  | -95.134+117.217i    -8.866-65.217i   -8.866+65.217i    -95.134-117.217i
   */
  {
    __attribute__ ((aligned (32))) complex double values[] = {
      209,
      -49,
      28-61*I,
      28+61*I,
     -95.13351365+117.21677477*I,
     0, //-8.86648635-65.21677477*I,
     0, //-8.86648635+65.21677477*I,
     0, //-95.13351365-117.21677477*I,
    };


    tft_inverse(values, 5, 3);

    assert_fp(values[0], 1*8);
    assert_fp(values[1], 34*8);
    assert_fp(values[2], 26*8);
    assert_fp(values[3], 95*8);
    assert_fp(values[4], 53*8);
    assert_fp(values[5], 0);
    assert_fp(values[6], 0);
    assert_fp(values[7], 0);
  }

  {
    __attribute__ ((aligned (32))) complex double values[] = {
        647.00000000+0*I,
        -29.00000000+0*I,
         81.00000000-74*I,
         81.00000000+74*I,
        -42.52186130+87.1543289*I,
        220.52186130-23.1543289*I,
        220.52186130+23.1543289*I,
        -42.52186130-87.1543289*I,
       -118.78486382+197.7215088*I,
        -48.89637704-63.4372376*I,
        -59.61723845+56.2080022*I,
         67.29847932+21.5077265*I,
       0, //  67.29847932-21.5077265*I,
       0, // -59.61723845-56.2080022*I,
       0, // -48.89637704+63.4372376*I,
       0, //-118.78486382-197.7215088*I,
    };


    tft_inverse(values, 12, 4);

    assert_fp(values[0],  51*16);
    assert_fp(values[1],  34*16);
    assert_fp(values[2],  26*16);
    assert_fp(values[3],  95*16);
    assert_fp(values[4],  53*16);
    assert_fp(values[5],  93*16);
    assert_fp(values[6],  41*16);
    assert_fp(values[7],  37*16);
    assert_fp(values[8],  91*16);
    assert_fp(values[9],   5*16);
    assert_fp(values[10], 47*16);
    assert_fp(values[11], 74*16);
    assert_fp(values[12],  0);
    assert_fp(values[13],  0);
    assert_fp(values[14],  0);
    assert_fp(values[15],  0);
  }
}
