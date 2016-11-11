#define _POSIX_C_SOURCE 199309L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // abort
#include <malloc.h>
#include <pmmintrin.h>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include <time.h>
#include "fft.h"
#include "bigfloat.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Helpers
double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

void dump_to_file(const char *path, const char* str, size_t len){
    FILE *file = fopen(path, "wb");

    if (file == NULL) {
      fprintf(stderr, "Cannot create file\n");
      abort();
    }

    fwrite(str, 1, len, file);
    fclose(file);
}

//  Returns the # of terms needed to reach a precision of p.
//
//  The taylor series converges to log(x!) / log(2) decimal digits after
//  x terms. So to find the number of terms needed to reach a precision of p
//  we need to solve this question for x:
//      p = log(x!) / log(4294967296)
//
//  This function solves this equation via binary search.
size_t e_terms(size_t p) {
  double sizeL = p * logl(4294967296) + 1;

  size_t a = 0;
  size_t b = 1;

  //  Double up
  while (lgammal(b) < sizeL) {
    b *= 2;
  }

  //  Binary search
  while (b - a > 1) {
    size_t m = (a + b) / 2;

    if (lgammal(m) < sizeL) {
      a = m;
    } else {
      b = m;
    }
  }

  return b + 2;
}

void e_BSR(bigfloat_t P, bigfloat_t Q, uint32_t a, uint32_t b) {
  //  Binary Splitting recursion for exp(1).

  if (b - a == 1){
    bigfloat_set(P, 1);
    bigfloat_set(Q, b);
    return;
  }

  uint32_t m = (a + b) / 2;

  bigfloat_t P0, Q0, P1, Q1;
  bigfloat_new(P0);
  bigfloat_new(Q0);
  bigfloat_new(P1);
  bigfloat_new(Q1);

  if (b - a < 1000) {
    //  No more threads.
    e_BSR(P0, Q0, a, m);
    e_BSR(P1, Q1, m, b);
  } else {
    cilk_spawn e_BSR(P0, Q0, a, m);
    e_BSR(P1, Q1, m, b);
    cilk_sync;
  }

  bigfloat_t tmp;
  bigfloat_new(tmp);
  bigfloat_mul(tmp, P0, Q1, 0);
  bigfloat_free(P0);
  bigfloat_add(P, tmp, P1, 0);
  bigfloat_free(tmp);
  bigfloat_free(P1);
  bigfloat_mul(Q, Q0, Q1, 0);
  bigfloat_free(Q1);
  bigfloat_free(Q0);
}

void e(size_t digits){
  //  The leading 2 doesn't count.
  digits++;

  size_t p = (digits + 8) / 9;
  size_t terms = e_terms(p);

  //  Limit Exceeded
  if ((uint32_t)terms != terms) {
    fprintf(stderr, "Limit Exceeded");
    abort();
  }

  double time0 = wall_clock();

  printf("Summing (%ld terms)...", terms);
  fflush(stdout);

  bigfloat_t P, Q;
  bigfloat_new(P);
  bigfloat_new(Q);
  e_BSR(P, Q, 0, (uint32_t)terms);
  double time1 = wall_clock();

  printf("ok [%f seconds]\n", time1 - time0);

  printf("Dividing...");
  fflush(stdout);

  bigfloat_t one;
  bigfloat_new(one);
  bigfloat_set(one, 1);
  bigfloat_t tmp;
  bigfloat_new(tmp);
  bigfloat_div(tmp, P, Q, p);
  bigfloat_free(Q);
  bigfloat_add(P, tmp, one, p);
  bigfloat_free(tmp);
  bigfloat_free(one);
  double time2 = wall_clock();
  printf("ok [%f seconds]\n", time2 - time1);

  printf("Writing hex digits...");
  fflush(stdout);
  size_t output_len = digits+2; // comma, first '2'
  char *out = (char*) malloc(output_len+1);
  size_t len = bigfloat_to_string(out, P, output_len, 16);
  dump_to_file("e-hex.txt", out, len);
  double time3 = wall_clock();
  printf("ok [%f seconds]\n", time3 - time2);

  printf("Converting to decimal...");
  fflush(stdout);
  bigfloat_t dec;
  bigfloat_new(dec);
  bigfloat_radix(dec, P);
  double time4 = wall_clock();
  printf("ok [%f seconds]\n", time4 - time3);

  printf("Writing decimal digits...");
  len = bigfloat_to_string(out, dec, output_len, 10);
  dump_to_file("e.txt", out, len);
  bigfloat_free(P);
  double time5 = wall_clock();
  printf("ok [%f seconds]\n", time5 - time4);
}

int main() {
  size_t digits = 1000;

  //  Determine minimum FFT size.
  size_t p      = 2*digits / 9 + 10;
  size_t k      = 0;
  size_t length = 1;

  while(length < 3*p) { // FIXME: Adjust to digits/point
    length *= 2;
    k++;
  }
  k++;

  printf("Configuration:\n");
  printf("Digits:           %ld\n", digits);
  printf("Threads:          %d\n", __cilkrts_get_nworkers());
  printf("CPU Features:     Using AVX\n");
  printf("Output File:      ./e.txt\n");
  printf("Max FFT required: 2^%ld\n", k);
  printf("\n");

  if(k>26) {
    k = 26;
  }
  double time1 = wall_clock();
  printf("Building FFT tables of size 2^%ld...", k);
  fflush(stdout);
  fft_ensure_table(k);
  double time2 = wall_clock();
  printf("ok [%f seconds]\n", time2 - time1);

  // Calculate e
  e(digits);
  double time3 = wall_clock();

  printf("\nCalculation complete [%f seconds]\n", time3 - time1);
}
