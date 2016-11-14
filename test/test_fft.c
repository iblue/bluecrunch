#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "fft.h"

#define assert_fp(a, b) _assert_fp(a, b, __FILE__, __LINE__)

static inline void _assert_fp(complex double a, complex double b, char* file, int line) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d\n",
      creal(b), cimag(b), creal(a), cimag(b), file, line);
    abort();
  }
}

int main() {
  fft_ensure_table(8); // up to 256 values

  {
    __attribute__ ((aligned (32))) complex double values[] = {1, 2, 3, 4};

    fft_forward(values, 4);

    assert_fp(values[0], 10);
    assert_fp(values[1], -2);
    assert_fp(values[2], -2-2*I);
    assert_fp(values[3], -2+2*I);
  }
}
