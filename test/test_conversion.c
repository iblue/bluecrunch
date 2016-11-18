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
  srand(17699823);

  for(size_t run=0;run<10000;run++) {
    uint32_t values[64];
    uint32_t target[64];
    complex double T[300];

    for(int i=0;i<64;i++) {
      values[i] = rand();
    }

    size_t len = rand()%64;
    size_t bpp = 7+rand()%(22-7);

    printf("Run %ld: Testing %ld values in %ld bit...", run, len, bpp);

    int_to_fft(T, 300, values, len, bpp);
    fft_to_int(T, 1, target, len, bpp); // 1 => do not scale


    for(int i=0;i<len;i++) {
      assert(values[i] == target[i]);
    }
    printf("[ok]\n");
  }
  return 0;
}
