#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // abort
#include <malloc.h>
#include <pmmintrin.h>

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "fft.h"
#include "bigfloat.h"
#include "bench.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Helpers
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

  task_start(1, "Summing");
  bigfloat_t P, Q;
  bigfloat_new(P);
  bigfloat_new(Q);
  e_BSR(P, Q, 0, (uint32_t)terms);
  task_end(1);


  task_start(1, "Dividing");
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
  task_end(1);

  task_start(1, "Writing hex digits");
  size_t output_len = digits+2; // comma, first '2'
  // We add 1 byte for the null terminator and 21 for overflows
  // We write 9 bytes at a time and stop when shit gets too large.
  // So just malloc a bit more and truncate, because I'm lazy.
  char *out = (char*) malloc(output_len+22);
  size_t len = bigfloat_to_string(out, P, output_len, 16);
  dump_to_file("e-hex.txt", out, len);
  task_end(1);

  task_start(1, "Converting to decimal");
  fflush(stdout);
  bigfloat_t dec;
  bigfloat_new(dec);
  bigfloat_radix(dec, P);
  task_end(1);

  task_start(1, "Writing decimal digits");
  len = bigfloat_to_string(out, dec, output_len, 10);
  dump_to_file("e.txt", out, len);
  bigfloat_free(P);
  task_end(1);
}

int main(int argc, char *argv[]) {
  if(argc != 2) {
    printf("Usage %s <number of digits> - Computes Eulers constant\n", argv[0]);
    exit(0);
  }

  size_t digits;
  sscanf(argv[1], "%ld", &digits);

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

  // Cap table sizes.
  if(k > 26) {
    k = 26;
  }

  task_start(0, "Starting calculation");
  task_start(1, "Precomputing twiddle factors");
  fft_ensure_table(3*length);
  task_end(1);

  // Calculate e
  e(digits);

  task_start(1, "Deallocating twiddle factors");
  fft_free_table();
  task_end(1);
  task_end(0);
}
