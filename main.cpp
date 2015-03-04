#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <pmmintrin.h>

#include <omp.h>

#include <string>

extern "C" {
  #include "fft.h"
}
#include "bigfloat.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Helpers
double wall_clock() {
  struct timespec t = {0,0};
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
//  The taylor series converges to log(x!) / log(10) decimal digits after
//  x terms. So to find the number of terms needed to reach a precision of p
//  we need to solve this question for x:
//      p = log(x!) / log(1000000000)
//
//  This function solves this equation via binary search.
size_t e_terms(size_t p) {
  double sizeL = p * logl(1000000000.0) + 1;

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

void e_BSR(BigFloat &P, BigFloat &Q, uint32_t a, uint32_t b, int tds = 1) {
  //  Binary Splitting recursion for exp(1).

  if (b - a == 1){
    bigfloat_set(P, 1, true);
    bigfloat_set(Q, b, true);
    return;
  }

  uint32_t m = (a + b) / 2;

  BigFloat P0, Q0, P1, Q1;
  bigfloat_new(P0);
  bigfloat_new(Q0);
  bigfloat_new(P1);
  bigfloat_new(Q1);

  if (b - a < 1000 || tds < 2) {
    //  No more threads.
    e_BSR(P0, Q0, a, m);
    e_BSR(P1, Q1, m, b);
  } else {
    //  Run sub-recursions in parallel.
    int tds0 = tds / 2;
    int tds1 = tds - tds0;
    #pragma omp parallel num_threads(2)
    {
      int tid = omp_get_thread_num();

      if (tid == 0){
        e_BSR(P0,Q0,a,m,tds0);
      }

      if (tid != 0 || omp_get_num_threads() < 2){
          e_BSR(P1,Q1,m,b,tds1);
      }
    }
  }

  BigFloat tmp;
  bigfloat_new(tmp);
  bigfloat_mul(tmp, P0, Q1, 0, tds);
  bigfloat_add(P, tmp, P1, 0);
  bigfloat_mul(Q, Q0, Q1, 0, tds);
  bigfloat_free(P0);
  bigfloat_free(P1);
  bigfloat_free(Q0);
  bigfloat_free(Q1);
  bigfloat_free(tmp);
}

void e(size_t digits, int tds){
  //  The leading 2 doesn't count.
  digits++;

  size_t p = (digits + 8) / 9;
  size_t terms = e_terms(p);

  //  Limit Exceeded
  if ((uint32_t)terms != terms) {
    fprintf(stderr, "Limit Exceeded");
    abort();
  }

  printf("Computing e...\n");
  printf("Algorithm: Taylor Series of exp(1)\n\n");

  double time0 = wall_clock();

  printf("Summing Series... %ld terms\n", terms);

  BigFloat P, Q;
  bigfloat_new(P);
  bigfloat_new(Q);
  e_BSR(P, Q, 0, (uint32_t)terms, tds);
  double time1 = wall_clock();

  printf("Time: %f\n", time1 - time0);

  printf("Division...\n");

  BigFloat one = BigFloat();
  bigfloat_new(one);
  bigfloat_set(one, 1, 1);
  BigFloat tmp;
  bigfloat_new(tmp);
  bigfloat_div(tmp, P, Q, p, tds);
  bigfloat_add(P, tmp, one, p);
  bigfloat_free(tmp);
  bigfloat_free(one);
  bigfloat_free(Q);
  double time2 = wall_clock();
  printf("Time: %f\n", time2 - time1);

  char *out;
  size_t len = bigfloat_to_string(out, P, digits);
  bigfloat_free(P);

  dump_to_file("e.txt", out, len);
}

int main() {
  // Enable threads and nested threading
  omp_set_nested(1);
  int threads = omp_get_max_threads();

  size_t digits = 10000000;

  //  Determine minimum FFT size.
  size_t p      = 2*digits / 9 + 10;
  size_t k      = 0;
  size_t length = 1;

  while(length < 5*p) {
    length *= 2;
    k++;
  }

  double time1 = wall_clock();
  printf("Bulding tables until size 2^%ld\n", k);

  // Precompute FFT twiddle factors
  fft_ensure_table(k);

  double time2 = wall_clock();
  printf("Time: %f\n", time2 - time1);

  // Calculate e
  e(digits, threads);
  double time3 = wall_clock();

  printf("Total Time = %f \n\n", time3 - time1);
}
