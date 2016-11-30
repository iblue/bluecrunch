/* gcc -Wall -Wextra -pedantic -march=native -mtune=native -O3 -std=c11 -o 01_basic 01_basic.c -lm */

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
void fft(complex double  *T, size_t len) {
  if (len == 1) {
    return;
  }

  size_t half_len = len / 2;

  double omega = 2 * M_PI / len;

  for (size_t c = 0; c < half_len; c++){
    double angle = omega * c;
    complex double twiddle_factor = cos(angle) + I*sin(angle);

    complex double a = T[c];
    complex double b = T[c + half_len];

    T[c           ] = a + b;
    T[c + half_len] = (a - b) * twiddle_factor;
  }

  fft(T + half_len, half_len);
  fft(T, half_len);
}

#define SIZE 1024
#define REP (100)

int main(void) {
  complex double *T = (complex double*)aligned_alloc(32, SIZE*sizeof(complex double));

  double start = wall_clock();
  for(int i=0;i<REP;i++) {
    fft(T, SIZE);
  }
  double end = wall_clock();
  printf("%f s\n", (end-start)/REP);
}
