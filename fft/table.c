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

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

my_complex_t* twiddle_table[32];
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
  my_complex_t *sub_table = (my_complex_t*)_mm_malloc(length*sizeof(my_complex_t), 32);

  for (size_t c = 0; c < length; c++){
      //  Generate Twiddle Factor
      double angle = omega * c;
      my_complex_t twiddle_factor;
      twiddle_factor.r = cos(angle);
      twiddle_factor.i = sin(angle);
      sub_table[c] = twiddle_factor;
  }

  //  Push into main table.
  twiddle_table[k] = sub_table;
  twiddle_table_size = k+1;
}

