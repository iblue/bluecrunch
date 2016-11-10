#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <pmmintrin.h>
#include <immintrin.h> // More Magic!
#include <omp.h>
#include "fft.h"
#include "intrinsic.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

void baileys_forward(complex double *T, int k, int tds) {
  size_t k1 = k / 2;
  size_t k2 = k - k1;
  size_t n1 = 1 << k1;
  size_t n2 = 1 << k2;
  size_t n  = 1 << k;

  // Do n1 transforms of size n2
  for(size_t i=0;i<n1;i++) {
    // FIXME: Parallelize
    fft_forward(T+i*n2, k1, tds);
  }

  // Do transpose and multiply by twiddle
  // FIXME: Not cache local. Can be done much faster
  for(size_t i=0;i<n1;i++) {
    for(size_t j=0;j<n2;j++) {
      // Swap
      complex double tmp = T[i+j*n1];
      T[i+j*n1] = T[j+i*n1];
      T[j+i*n1] = tmp;
    }
  }

  double omega = 2 * M_PI / n;

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
    fft_forward(T+i*n1, k2, tds);
  }
}
