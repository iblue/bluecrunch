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
#include "intrinsic.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

void fft_forward(complex double *T, int k){
  // (Bit reversed) 2-point DFT
  if(k == 1) {
    dft_2p(T);
    return;
  }

  size_t length      = 1 << k;
  size_t half_length = length / 2;

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[k];

  for (size_t n = 0; n < half_length; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddles = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    fft_forward_butterfly(twiddles, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }

  if(k >= FFT_THRESHOLD_K) {
    cilk_spawn fft_forward(T, k - 1);
    fft_forward(T + half_length, k - 1);
    cilk_sync;
  } else {
    fft_forward(T, k - 1);
    fft_forward(T + half_length, k - 1);
  }
}

// Runs FFT without cached twiddles
// FIXME: Merge?
void fft_forward_uncached(complex double *T, int k){
  //  End recursion at 2 points.
  if (k == 1){
    dft_2p(T);
    return;
  }

  size_t length = 1 << k;
  double omega = 2 * M_PI / length;
  size_t half_length = length / 2;

  //  Perform FFT reduction into two halves.
  for (size_t c = 0; c < half_length; c++){
    //  Grab Twiddle Factor
    double angle = omega*c;
    double r = cos(angle);
    double i = sin(angle);
    __m128d r0 = _mm_loaddup_pd(&r);
    __m128d i0 = _mm_loaddup_pd(&i);

    //  Grab elements
    __m128d a0 = ((__m128d*)T)[c];
    __m128d b0 = ((__m128d*)T)[c + half_length];

    //  Perform butterfly
    __m128d c0,d0;
    c0 = _mm_add_pd(a0,b0);
    d0 = _mm_sub_pd(a0,b0);

    ((__m128d*)T)[c] = c0;

    //  Multiply by twiddle factor.
    c0 = _mm_mul_pd(d0,r0);
    d0 = _mm_mul_pd(_mm_shuffle_pd(d0,d0,1),i0);
    c0 = _mm_addsub_pd(c0,d0);

    ((__m128d*)T)[c + half_length] = c0;
  }

  if(k-1 < twiddle_table_size) {
    if(k >= FFT_THRESHOLD_K) {
      cilk_spawn fft_forward(T, k - 1);
      fft_forward(T + half_length, k - 1);
      cilk_sync;
    } else {
      fft_forward(T, k - 1);
      fft_forward(T + half_length, k - 1);
    }
  } else {
    if(k >= FFT_THRESHOLD_K) {
      cilk_spawn fft_forward_uncached(T, k - 1);
      fft_forward_uncached(T + half_length, k - 1);
      cilk_sync;
    } else {
      fft_forward_uncached(T, k - 1);
      fft_forward_uncached(T + half_length, k - 1);
    }
  }
}

void fft_inverse(complex double *T, int k){
  //  Fast Fourier Transform
  //  This function performs an inverse FFT of length 2^k.

  //  This is a Decimation-in-Time (DIT) FFT.
  //  The frequency domain input must be in bit-reversed order.

  //Parameters:
  //  -   T           -   Pointer to array.
  //  -   k           -   2^k is the size of the transform

  //  End recursion at 2 points.
  if (k == 1){
    dft_2p(T);
    return;
  }

  size_t length = 1 << k;
  size_t half_length = length / 2;

  if(k >= FFT_THRESHOLD_K) {
    cilk_spawn fft_inverse(T, k - 1);
    fft_inverse(T + half_length, k - 1);
    cilk_sync;
  } else {
    fft_inverse(T, k - 1);
    fft_inverse(T + half_length, k - 1);
  }

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[k];

  //  Perform FFT reduction into two halves.
  for (size_t n = 0; n < half_length; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddle = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    fft_inverse_butterfly(twiddle, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }
}

void fft_inverse_uncached(complex double *T, int k){
  //  Fast Fourier Transform
  //  This function performs an inverse FFT of length 2^k.

  //  This is a Decimation-in-Time (DIT) FFT.
  //  The frequency domain input must be in bit-reversed order.

  //Parameters:
  //  -   T           -   Pointer to array.
  //  -   k           -   2^k is the size of the transform

  //  End recursion at 2 points.
  if (k == 1){
    dft_2p(T);
    return;
  }

  size_t length = 1 << k;
  double omega = 2 * M_PI / length;
  size_t half_length = length / 2;

  if(k - 1 < twiddle_table_size) {
    if(k >= FFT_THRESHOLD_K) {
      cilk_spawn fft_inverse(T, k - 1);
      fft_inverse(T + half_length, k - 1);
      cilk_sync;
    } else {
      fft_inverse(T, k - 1);
      fft_inverse(T + half_length, k - 1);
    }
  } else {
    if(k >= FFT_THRESHOLD_K) {
      cilk_spawn fft_inverse_uncached(T, k - 1);
      fft_inverse_uncached(T + half_length, k - 1);
      cilk_sync;
    } else {
      fft_inverse_uncached(T, k - 1);
      fft_inverse_uncached(T + half_length, k - 1);
    }
  }

  //  Perform FFT reduction into two halves.
  for (size_t c = 0; c < half_length; c++){
    //  Generate Twiddle Factor
    double angle = omega * c;
    double r = cos(angle);
    double i = sin(angle);
    __m128d r0 = _mm_loaddup_pd(&r);
    __m128d i0 = _mm_loaddup_pd(&i);
    i0 = _mm_xor_pd(i0,_mm_set1_pd(-0.0));

    //  Grab elements
    __m128d a0 = ((__m128d*)T)[c];
    __m128d b0 = ((__m128d*)T)[c + half_length];

    //  Perform butterfly
    __m128d c0,d0;

    //  Multiply by twiddle factor.
    c0 = _mm_mul_pd(b0,r0);
    d0 = _mm_mul_pd(_mm_shuffle_pd(b0,b0,1),i0);
    c0 = _mm_addsub_pd(c0,d0);

    b0 = _mm_add_pd(a0,c0);
    d0 = _mm_sub_pd(a0,c0);

    ((__m128d*)T)[c] = b0;
    ((__m128d*)T)[c + half_length] = d0;
  }
}

