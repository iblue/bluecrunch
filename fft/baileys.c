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

static inline void swap(complex double* a, complex double* b) {
  complex double tmp = *a;
  *a = *b;
  *b = tmp;
}

static inline void transpose(complex double *m, size_t n) {
  size_t block = 0;
  size_t size  = 8;

  for(block = 0; block + size - 1 < n; block += size){
    for(size_t i = block; i < block+size; i++) {
      for(size_t j = i + 1; j < block + size; j++) {
        swap(&m[i*n + j], &m[j*n + i]);
      }
    }
    for(size_t i = block + size; i < n; i++) {
      for(size_t j = block; j < block + size; j++) {
        swap(&m[i*n + j], &m[j*n + i]);
      }
    }
  }

  for(size_t i=block;i<n;++i) {
    for(size_t j=i+1;j<n;++j) {
      swap(&m[i*n + j], &m[j*n + i]);
    }
  }
}

void baileys_forward(complex double *T, size_t length) {
  // FFT conform sqrt for dummies :)
  size_t n1 = 1;
  size_t n2 = length;

  while(n1 < n2) {
    n2 /= 2;
    n1 *= 2;
  }

  transpose(T, n1);

  // Do n1 transforms of size n2
  for(size_t i=0;i<n1;i++) {
    // FIXME: Parallelize
    fft_forward(T+i*n2, n2);
  }

  transpose(T, n2);

  double omega = 2*M_PI / length;

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
}
