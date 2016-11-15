#include <stdio.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include "fft.h"

#define assert_fp(a, b) _assert_fp(a, b, __FILE__, __LINE__)
static inline void _assert_fp(complex double a, complex double b, char* file, int line) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d\n",
      creal(b), cimag(b), creal(a), cimag(a), file, line);
    abort();
  }
}

int main(void) {
  // 12 bits, 1/3 words
  {
    complex double T[16];

    uint32_t words[] = {0x1234abcd};

    int_to_fft(T, 4, words, 1, 12);

    assert_fp(T[0], 0xbcd);
    assert_fp(T[1], 0x34a);
    assert_fp(T[2], 0x12);

    for(size_t i=3;i<16;i++) {
      assert_fp(T[i], 0x0);
    }
  }

  // 12 bits, 2/3 words
  {
    complex double T[16];

    uint32_t words[] = {0x1234abcd, 0xdeadbeef};

    int_to_fft(T, 4, words, 2, 12);

    assert_fp(T[0], 0xbcd);
    assert_fp(T[1], 0x34a);
    assert_fp(T[2], 0xf12);
    assert_fp(T[3], 0xbee);
    assert_fp(T[4], 0xead);
    assert_fp(T[5], 0xd);

    for(size_t i=6;i<16;i++) {
      assert_fp(T[i], 0x0);
    }
  }

  // 12 bits, 3/3 words
  {
    complex double T[16];

    uint32_t words[] = {0x1234abcd, 0xdeadbeef, 0xcafebabe};

    int_to_fft(T, 4, words, 3, 12);

    assert_fp(T[0], 0xbcd);
    assert_fp(T[1], 0x34a);
    assert_fp(T[2], 0xf12);
    assert_fp(T[3], 0xbee);
    assert_fp(T[4], 0xead);
    assert_fp(T[5], 0xbed);
    assert_fp(T[6], 0xeba);
    assert_fp(T[7], 0xcaf);

    for(size_t i=8;i<16;i++) {
      assert_fp(T[i], 0x0);
    }
  }

  // 8 bits, 1/1 words
  {
    complex double T[16];

    uint32_t words[] = {0x1234abcd};

    int_to_fft(T, 4, words, 1, 8);

    assert_fp(T[0], 0xcd);
    assert_fp(T[1], 0xab);
    assert_fp(T[2], 0x34);
    assert_fp(T[3], 0x12);

    for(size_t i=4;i<16;i++) {
      assert_fp(T[i], 0x0);
    }
  }

  return 0;
}
