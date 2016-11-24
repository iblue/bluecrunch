#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // abort
#include <malloc.h>
#include <string.h>
#include <assert.h>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "fft.h"
#include "bench.h"

#define assert_fp(a, b, idx) _assert_fp(a, b, idx, __FILE__, __LINE__)
static inline void _assert_fp(complex double a, complex double b, size_t idx, char* file, int line) {
  double radius = 1e-10;

  if(cabs(a-b) > radius) {
    fprintf(stderr, "Expected: (%f + %fi) but got (%f + %fi) in %s:%d (iteration %ld)\n",
      creal(b), cimag(b), creal(a), cimag(a), file, line, idx);
    abort();
  }
}

#define K (20)

int main(int argc, char *argv[]) {
  // Test FFT performance up to 2^26
  complex double *T    = (complex double*)_mm_malloc((1 << K)*sizeof(complex double), 32);
  complex double *vals = (complex double*)_mm_malloc((1 << K)*sizeof(complex double), 32);

  {
    printf("Initializing values...");
    fflush(stdout);

    double start = wall_clock();
    for(size_t i=0;i<(1<<K);i++) {
      //T[i] = ((double)rand()/(double)RAND_MAX) + I*((double)rand()/(double)RAND_MAX);
      T[i] = 0.0;
    }

    memcpy(vals, T, (1<<K)*sizeof(complex double));

    for(size_t i=0;i<(1<<K);i++) {
      assert_fp(T[i], vals[i], i);
    }
    double end = wall_clock();
    printf(" -> %f seconds\n", end - start);
  }

  {
    printf("Generating tables...");
    fflush(stdout);

    double start = wall_clock();
    fft_ensure_table(1 << K);
    double end = wall_clock();
    printf(" -> %f seconds\n", end - start);
  }

  printf("== Old algorithm ==\n");
  for(size_t i=1;i<=K;i++) {
    printf("Benchmarking FFT of 2^%ld", i);
    fflush(stdout);
    double start = wall_clock();
    fft_forward(T, 1 << i);
    fft_inverse(T, 1 << i);
    double end = wall_clock();
    /*
    for(size_t i=1;i<=(1<<i);i++) {
      assert_fp(T[i], vals[i], i);
    }
    */

    printf(" -> %f seconds\n", end - start);
  }

  printf("== New algorithmn ==\n");
  for(size_t i=1;i<=K;i++) {
    printf("Benchmarking FFT of 2^%ld", i);
    fflush(stdout);
    double start = wall_clock();
    fftnew_forward(T, 1 << i);
    fftnew_inverse(T, 1 << i);
    double end = wall_clock();
    for(size_t i=1;i<=(1<<i);i++) {
      assert_fp(T[i], vals[i], i);
    }
    printf(" -> %f seconds\n", end - start);
  }
}
