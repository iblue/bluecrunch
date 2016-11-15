#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <pmmintrin.h>
#include <immintrin.h> // More Magic!
#include "fft.h"
#include "intrinsic.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

static inline void transpose(complex double *T, size_t n1, size_t n2) {
  // FIXME: Not cache local. Can be done much faster
  // FIXME: Not memory efficient
  complex double *N = malloc(sizeof(complex double)*(n1*n2));

  for(size_t i=0;i<n1;i++) {
    for(size_t j=0;j<n2;j++) {
      N[i+j*n1] = T[j+i*n2];
    }
  }

  memcpy(T, N, sizeof(complex double)*(n1*n2));
  free(N);

}

void baileys_forward(complex double *T, size_t length) {
  // FFT conform sqrt for dummies :)
  size_t n1 = 1;
  size_t n2 = length;
  while(n1 < n2) {
    n2 /= 2;
    n1 *= 2;
  }

  transpose(T, n1, n2);

  // Do n1 transforms of size n2
  for(size_t i=0;i<n1;i++) {
    // FIXME: Parallelize
    fft_forward(T+i*n2, n2);
  }

  double omega = 2 * M_PI / length;

  transpose(T, n2, n1);

  // FIXME: Do transpose and multiply in one step
  // FIXME: Use precomputed twiddles if available
  for(size_t i=0;i<n1;i++) {
    for(size_t j=0;j<n2;j++) {
      double angle = omega * (i*j);
      complex double twiddle_factor = cos(angle)+I*sin(angle);

      T[i+j*n1] *= twiddle_factor;
    }
  }

  // Do n2 transforms of size n1
  for(size_t i=0;i<n2;i++) {
    // FIXME: Parallelize
    fft_forward(T+i*n1, n1);
  }

  transpose(T, n1, n2);
}
