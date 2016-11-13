#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <pmmintrin.h>
#include <immintrin.h> // More Magic!
#include <cilk/cilk.h>
#include "fft.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

complex double* twiddle_table[32];
int twiddle_table_size = 0;

void fft_ensure_table(int k) {
  //  Makes sure the twiddle factor table is large enough to handle an FFT of
  //  size 2^k.

  //  Do one level at a time
  if (k - 1 > 0) {
    fft_ensure_table(k - 1);
  }

  size_t length = 1 << k;
  double omega = 2 * M_PI / length;
  length /= 2;

  //  Build the sub-table.
  complex double *sub_table = (complex double*)_mm_malloc(length*sizeof(complex double), 32);

  cilk_for(size_t c = 0; c < length; c++) {
    //  Generate Twiddle Factor
    double angle = omega * c;
    complex double twiddle_factor = cos(angle)+I*sin(angle);
    sub_table[c] = twiddle_factor;
  }

  //  Push into main table.
  twiddle_table[k] = sub_table;
  twiddle_table_size = k+1;
}

void fft_free_table() {
  for(size_t i=0;i<twiddle_table_size;i++) {
    _mm_free(twiddle_table[i]);
  }
}
