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
#include <x86intrin.h>

double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

#define SIZE (1024*1024*64)
#define REP (100)

void __attribute__ ((noinline)) seq(double *mem) {
  for(size_t i=0;i<SIZE;i++) {
    mem[i] *= 3.7f;
  }
}

void __attribute__ ((noinline)) opt(double *mem) {
  static const double mul[4] = {3.7f, 3.7f, 3.7f, 3.7f};

  __m256d xmm0 = _mm256_load_pd(mul);

  for(size_t i=0;i<SIZE;i+=32) {
    __m256d xmm1 = _mm256_load_pd(mem+i);
    __m256d xmm2 = _mm256_load_pd(mem+i+4);
    __m256d xmm3 = _mm256_load_pd(mem+i+8);
    __m256d xmm4 = _mm256_load_pd(mem+i+12);
    __m256d xmm5 = _mm256_load_pd(mem+i+16);
    __m256d xmm6 = _mm256_load_pd(mem+i+20);
    __m256d xmm7 = _mm256_load_pd(mem+i+24);
    __m256d xmm8 = _mm256_load_pd(mem+i+28);
    __m256d xmm9 = _mm256_load_pd(mem+i+32);

    xmm1 = _mm256_mul_pd(xmm0, xmm1);
    xmm2 = _mm256_mul_pd(xmm0, xmm2);
    xmm3 = _mm256_mul_pd(xmm0, xmm3);
    xmm4 = _mm256_mul_pd(xmm0, xmm4);
    xmm5 = _mm256_mul_pd(xmm0, xmm5);
    xmm6 = _mm256_mul_pd(xmm0, xmm6);
    xmm7 = _mm256_mul_pd(xmm0, xmm7);
    xmm8 = _mm256_mul_pd(xmm0, xmm8);
    xmm9 = _mm256_mul_pd(xmm0, xmm9);

    _mm256_store_pd(mem+i,    xmm1);
    _mm256_store_pd(mem+i+4,  xmm2);
    _mm256_store_pd(mem+i+8,  xmm3);
    _mm256_store_pd(mem+i+12, xmm4);
    _mm256_store_pd(mem+i+16, xmm5);
    _mm256_store_pd(mem+i+20, xmm6);
    _mm256_store_pd(mem+i+24, xmm7);
    _mm256_store_pd(mem+i+28, xmm8);
    _mm256_store_pd(mem+i+32, xmm9);
  }
}

void __attribute__ ((noinline)) opt2(double *mem) {
  //__builtin_prefetch(mem);
  //__builtin_prefetch(mem+32);
  static const double mul[4] = {3.7f, 3.7f, 3.7f, 3.7f};

  __m256d xmm0 = _mm256_load_pd(mul);

  for(size_t i=0;i<SIZE;i+=32) {
    __m256d xmm1 = _mm256_load_pd(mem+i);
    xmm1 = _mm256_mul_pd(xmm0, xmm1);
    _mm256_store_pd(mem+i,    xmm1);

    __m256d xmm2 = _mm256_load_pd(mem+i+4);
    xmm2 = _mm256_mul_pd(xmm0, xmm2);
    _mm256_store_pd(mem+i+4,  xmm2);

    __m256d xmm3 = _mm256_load_pd(mem+i+8);
    xmm3 = _mm256_mul_pd(xmm0, xmm3);
    _mm256_store_pd(mem+i+8,  xmm3);

    __m256d xmm4 = _mm256_load_pd(mem+i+12);
    xmm4 = _mm256_mul_pd(xmm0, xmm4);
    _mm256_store_pd(mem+i+12, xmm4);

    __m256d xmm5 = _mm256_load_pd(mem+i+16);
    xmm5 = _mm256_mul_pd(xmm0, xmm5);
    _mm256_store_pd(mem+i+16, xmm5);

    __m256d xmm6 = _mm256_load_pd(mem+i+20);
    xmm6 = _mm256_mul_pd(xmm0, xmm6);
    _mm256_store_pd(mem+i+20, xmm6);

    __m256d xmm7 = _mm256_load_pd(mem+i+24);
    xmm7 = _mm256_mul_pd(xmm0, xmm7);
    _mm256_store_pd(mem+i+24, xmm7);

    __m256d xmm8 = _mm256_load_pd(mem+i+28);
    xmm8 = _mm256_mul_pd(xmm0, xmm8);
    _mm256_store_pd(mem+i+28, xmm8);

    __m256d xmm9 = _mm256_load_pd(mem+i+32);
    xmm9 = _mm256_mul_pd(xmm0, xmm9);
    _mm256_store_pd(mem+i+32, xmm9);
  }
}

void __attribute__ ((noinline)) opt2_sse(double *mem) {
  __builtin_prefetch(mem);
  __builtin_prefetch(mem, 1);
  //__builtin_prefetch(mem);
  //__builtin_prefetch(mem+32);
  static const double mul[4] = {3.7f, 3.7f};

  __m128d xmm0 = _mm_load_pd(mul);

  for(size_t i=0;i<SIZE;i+=32) {
    __builtin_prefetch(mem+i+32);
    __builtin_prefetch(mem+i+32, 1);
    __m128d xmm1 = _mm_load_pd(mem+i);
    xmm1 = _mm_mul_pd(xmm0, xmm1);
    _mm_store_pd(mem+i,    xmm1);

    __m128d xmm2 = _mm_load_pd(mem+i+4);
    xmm2 = _mm_mul_pd(xmm0, xmm2);
    _mm_store_pd(mem+i+4,  xmm2);

    __m128d xmm3 = _mm_load_pd(mem+i+8);
    xmm3 = _mm_mul_pd(xmm0, xmm3);
    _mm_store_pd(mem+i+8,  xmm3);

    __m128d xmm4 = _mm_load_pd(mem+i+12);
    xmm4 = _mm_mul_pd(xmm0, xmm4);
    _mm_store_pd(mem+i+12, xmm4);

    __m128d xmm5 = _mm_load_pd(mem+i+16);
    xmm5 = _mm_mul_pd(xmm0, xmm5);
    _mm_store_pd(mem+i+16, xmm5);

    __m128d xmm6 = _mm_load_pd(mem+i+20);
    xmm6 = _mm_mul_pd(xmm0, xmm6);
    _mm_store_pd(mem+i+20, xmm6);

    __m128d xmm7 = _mm_load_pd(mem+i+24);
    xmm7 = _mm_mul_pd(xmm0, xmm7);
    _mm_store_pd(mem+i+24, xmm7);

    __m128d xmm8 = _mm_load_pd(mem+i+28);
    xmm8 = _mm_mul_pd(xmm0, xmm8);
    _mm_store_pd(mem+i+28, xmm8);

    __m128d xmm9 = _mm_load_pd(mem+i+32);
    xmm9 = _mm_mul_pd(xmm0, xmm9);
    _mm_store_pd(mem+i+32, xmm9);
  }
}

void linewise(double *mem) {
  for(size_t i=0;i<SIZE;i+=8) {
    mem[i]   *= 3.7f;
    mem[i+1] *= 3.7f;
    mem[i+2] *= 3.7f;
    mem[i+3] *= 3.7f;
    mem[i+4] *= 3.7f;
    mem[i+5] *= 3.7f;
    mem[i+6] *= 3.7f;
    mem[i+7] *= 3.7f;
  }
}

void interleaved(double *mem) {
  for(size_t i=0;i<SIZE;i+=16) {
    mem[i]   *= 3.7f;
    mem[i+1] *= 3.7f;
    mem[i+2] *= 3.7f;
    mem[i+3] *= 3.7f;
    mem[i+4] *= 3.7f;
    mem[i+5] *= 3.7f;
    mem[i+6] *= 3.7f;
    mem[i+7] *= 3.7f;
  }
  for(size_t i=0;i<SIZE;i+=16) {
    mem[i+8]  *= 3.7f;
    mem[i+9]  *= 3.7f;
    mem[i+10] *= 3.7f;
    mem[i+11] *= 3.7f;
    mem[i+12] *= 3.7f;
    mem[i+13] *= 3.7f;
    mem[i+14] *= 3.7f;
    mem[i+15] *= 3.7f;
  }
}

int main(void) {
  double start,end;
  double *mem = (double*)aligned_alloc(64, SIZE*sizeof(double));
  memset(mem, 0, SIZE*sizeof(double));

  // Access sequential
  start = wall_clock();
  seq(mem);
  end = wall_clock();
  printf("sequential: %f s\n", end-start);

  // Access sequential
  start = wall_clock();
  opt(mem);
  end = wall_clock();
  printf("avx: %f s\n", end-start);

  // Access sequential
  start = wall_clock();
  opt2(mem);
  end = wall_clock();
  printf("opt2: %f s\n", end-start);

  // Access sequential
  start = wall_clock();
  opt2_sse(mem);
  end = wall_clock();
  printf("opt2_sse: %f s\n", end-start);

  // Access strided (64 byte cache line);
  start = wall_clock();
  linewise(mem);
  end = wall_clock();
  printf("linewise: %f s\n", end-start);

  // Access interleaved
  start = wall_clock();
  interleaved(mem);
  end = wall_clock();
  printf("interleaved: %f s\n", end-start);
}
