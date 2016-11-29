#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 700
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <complex.h>

double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

// Amount of doubles that fits into L2 Cache
// In L1, we can fit 256 rows * 64 byte cache line (and we may need padding of
// 64 Bytes per row to break associativity).
// 256 rows is the maximum!
// There is no limit for the cols, but if we sub-divide, we may need additional padding.
#define ROWS (256)
#define COLS (512*1024)
#define SIZE (ROWS*COLS) // 256 rows x 512 cols
#define REP (1)

// 1 complex double = 16 bytes
// Bulldozer:
// L1: 16K 4-way associative with 64 Byte Cache lines (per Core) (256 lines, each containing 4 doubles)
// L2:  2M 16-way associative with 64 Byte Cache lines (shared among 2 Cores) (32768 lines)
// L3:  8M 64-way associative with 64 Byte Cache lines (shared among all Cores) (131072 lines)
// 1 Cache Line = 64 Bytes = 8 doubles = 4 complex doubles
// To use all lines, we need 64 B Padding in each row.
// AVX registers: 512 bytes = 8 Cache Lines

// Simulates normal fft access pattern
complex double fft(complex double *T, size_t len) {
  complex double sum = 0.0;

  if(len == 2) {
    sum = T[0] + T[1];
    return sum;
  }

  for(size_t i=0;i<len;i++) {
    sum += T[i];
  }

  sum += fft(T, len/2);
  sum += fft(T+len/2, len/2);

  return sum;
}

complex double strided_fft(complex double *T, size_t col, size_t rows) {
  if(rows == 1) {
    return T[0] + T[1] + T[2] + T[3];
  }

  complex double sum = 0.0;
  // FIXME: More realistic access pattern
  //
  // 0 1 2 3
  // 512 513 514 515
  // ...
  // 130560
  for(size_t row=0;row<rows;row++) {
    sum += T[col+row*COLS] + T[col+row*COLS+1] + T[col+row*COLS+2] + T[col+row*COLS+3];
  }

  sum += strided_fft(T, col, rows/2);
  sum += strided_fft(T+(rows/2)*COLS, col, rows/2);

  return sum;
}

// Simulates access patterns for 1M FFT
int main(void) {
  double start, end;

  // Align to cache lines
  complex double *T = aligned_alloc(64, SIZE*sizeof(complex double));
  memset(T, 0, SIZE*sizeof(complex double));

  // 1. Do FFT on the cols (vectorized)
  // 256 rows => 512 cols
  //
  // row 0: 0..511 ( read 4 at a time)
  // row 1: 512..1023
  // ...
  // row 255: 130560..131072

  complex double sum = 0.0f;

  start = wall_clock();
  for(size_t col=0;col<COLS;col+=4) {
    sum += strided_fft(T, col, ROWS);
  }

  // Now do a sequential read
  for(size_t row=0;row<ROWS;row++) {
    sum += fft(T+COLS*row, COLS);
  }
  end = wall_clock();

  printf("%f + %f*I\n", creal(sum), cimag(sum));

  printf("Improved access pattern: %f s\n", end - start);


  start = wall_clock();
  sum = fft(T, SIZE);
  end = wall_clock();
  printf("%f + %f*I\n", creal(sum), cimag(sum));
  printf("Normal access pattern: %f s\n", end - start);
}
