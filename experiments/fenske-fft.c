#include <complex.h> // complex
#include <stddef.h>  // size_t
#include <stdio.h>   // printf
#include <stdlib.h>  // rand, aligned_alloc

#define L1D_CACHE_LINE_SIZE (64)      // L1 Cache line size in bytes (64 B)
#define L1D_CACHE_SIZE      (16384)   // L1 Total Cache size in bytes (16 KB)

#define L2D_CACHE_LINE_SIZE (1024)    // L2 Cache line size in bytes (1 KB)
#define L2D_CACHE_SIZE      (2097152) // L2 Total Cache size in bytes (2 MB)

#define L3D_CACHE_LINE_SIZE (1024)    // L3 Cache line size in bytes (1 KB)
#define L3D_CACHE_SIZE      (8388608) // L3 Total Cache size in bytes (8 MB)

static inline void swap(complex double* a, complex double* b) {
  complex double tmp = *a;
  *a = *b;
  *b = tmp;
}

void transpose(complex double *m, size_t n) {
  size_t block = 0;
  size_t size  = 8;

  for(block = 0; block + size - 1 < n; block += size){
    for(size_t i = block; i < block+size; i++) {
      for(size_t j = i + 1; j < block + size; j++) {
        swap(&m[i*n + j], &m[j*n + i]);
      }
    }
    for(size_t i = block + size; i < n; i++) {
      for(size_t j = block; j < block + size; j++) {
        swap(&m[i*n + j], &m[j*n + i]);
      }
    }
  }

  for(size_t i=block;i<n;++i) {
    for(size_t j=i+1;j<n;++j) {
      swap(&m[i*n + j], &m[j*n + i]);
    }
  }
}

void fenske_fft(complex double *values, size_t len) {
  size_t n1 = 1;
  size_t n2 = len;

  while(n1 < n2) {
    n1 *= 2;
    n2 /= 2;
  }

  transpose(values, n1);

  // FFT on rows
  for(size_t i = 0; i < n1; i++) {
    fft_forward(values+i*n2, n2);
  }

  // Multiply by twiddles

  transpose(values, n2);

  // FFT on rows
  for(size_t i = 0; i < n2; i++) {
    fft_forward(values+i*n1, n1);
  }

}

#define FFT_SIZE 1024

int main(void) {
  complex double *values = aligned_alloc(64, FFT_SIZE*sizeof(complex double)); // align on cache lines.

  for(size_t i=0;i<FFT_SIZE;i++) {
    values[i] = (double)rand()/(double)(RAND_MAX);
  }

  fenske_fft(values, FFT_SIZE);

  return 0;
}
