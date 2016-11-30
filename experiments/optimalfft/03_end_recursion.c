/* gcc -Wall -Wextra -pedantic -march=native -mtune=native -O3 -std=c11 -o 03_end_recursion 03_end_recursion.c -lm */

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
#include <math.h>

double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define K (10)
#define SIZE (1<<K)
#define REP (10000)

complex double *table[K];

void fft(complex double  *T, size_t len) {
  // Gives 18% speedup over recursion end at 1
  if (len == 4) {
    complex double a = T[0];
    complex double b = T[2];
    complex double c = T[1];
    complex double d = T[3];
    T[0] = a+b+c+d;
    T[1] = a-b+c-d;
    T[2] = a+b*I-c-d*I;
    T[3] = a-b*I-c+d*I;
    return;
  }

  size_t half_len = len / 2;

  size_t tbl_idx = K;
  while(1) {
    if(len == 1) {
      break;
    }
    len /= 2;
    tbl_idx--;
  }

  complex double *tbl = table[tbl_idx];

  for (size_t c = 0; c < half_len; c++){
    complex double twiddle_factor = tbl[c];

    complex double a = T[c];
    complex double b = T[c + half_len];

    T[c           ] = a + b;
    T[c + half_len] = (a - b) * twiddle_factor;
  }

  fft(T + half_len, half_len);
  fft(T, half_len);
}

void precompute(complex double **table, size_t len) {
  size_t c=0;
  while(1) {
    len /= 2;
    double omega = 2 * M_PI / len;
    table[c] = (complex double*)aligned_alloc(32, len*sizeof(complex double));
    for(size_t i=0;i<len;i++) {
      double angle = omega * i;
      complex double twiddle_factor = cos(angle) + I*sin(angle);
      table[c][i] = twiddle_factor;
    }
    if(len == 1) {
      break;
    }
    c++;
  }
}

int main(void) {
  complex double *T = (complex double*)aligned_alloc(32, SIZE*sizeof(complex double));

  double start = wall_clock();
  precompute(table, SIZE);
  for(int i=0;i<REP;i++) {
    fft(T, SIZE);
  }
  double end = wall_clock();
  printf("%f us\n", (end-start)*1000000/REP);
}
