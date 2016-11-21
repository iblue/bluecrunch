#include <stdio.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include "fft.h"
#include "bench.h"

#define assert_fp(a, b) _assert_fp(a, b, __FILE__, __LINE__)
static inline void _assert_fp(complex double a, complex double b, char* file, int line) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d\n",
      creal(b), cimag(b), creal(a), cimag(a), file, line);
    abort();
  }
}

#define K (26)

int main(void) {
  // Test FFT performance up to 2^26
  complex double *T = (complex double*)_mm_malloc((1 << K)*sizeof(complex double), 32);

  fft_ensure_table(1 << K);

  for(size_t i=1;i<=K;i++) {
    printf("Benchmarking FFT of 2^%ld", i);
    double start = wall_clock();
    fft_forward(T, 1 << i);
    double end = wall_clock();
    printf(" -> %f seconds\n", end - start);
  }
}
