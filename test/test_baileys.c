#include "fft.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define assert_fp(a, b) _assert_fp(a, b, __FILE__, __LINE__)
static inline void _assert_fp(complex double a, complex double b, char* file, int line) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d\n",
      creal(b), cimag(b), creal(a), cimag(a), file, line);
    abort();
  }
}

int main() {
  fft_ensure_table(8); // up to 256 values

  {
    __attribute__ ((aligned (32))) complex double values1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    __attribute__ ((aligned (32))) complex double values2[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    baileys_forward(values1, 16);
    fft_forward(values2, 16);

    // FIXME: FP accuracy!
    for(size_t i=0;i<16;i++) {
      assert_fp(values1[i], values2[i]);
    }
  }
}
