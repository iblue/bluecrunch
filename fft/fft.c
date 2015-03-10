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

// Runs a forward butterfly on 2 elements simultaniously
// Arguments:
// - twiddles: Vektor of twiddles [r0, i0, r1, i1]
// - T0:       Pointer to 2 complex doubles (first half)
// - T1:       Pointer to 2 complex doubles (second half)
//
//  Perform FFT reduction into two halves.
// Input: a[0], a[1], b[0], b[1], W[0], W[1]
//        ----------  ----------  ----------
//            |           |           |
//      real & imag T[c]  |    real & imag root
//               real & imag T[c+N/2]
//
// Output:
//  a[0] <- a[0] + b[0]
//  a[1] <- a[1] + b[1]
//  b[0] <- W[0]*(a[0] - b[0]) - W[1]*(a[1] - b[1])
//  b[1] <- W[0]*(a[1] - b[1]) - W[1]*(a[0] - b[0])
static inline void fft_forward_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
  __m256d r   = _mm256_unpacklo_pd(twiddles, twiddles);  //   r = [r0,r0,r1,r1]
  __m256d i   = _mm256_unpackhi_pd(twiddles, twiddles);  //   i = [i0,i0,i1,i1]

  //  Grab elements
  __m256d a = _mm256_load_pd((double*)T0);  // a = [a0r, a0i, a1r, a1i]
  __m256d b = _mm256_load_pd((double*)T1);  // b = [b0r, b0i, b1r, b1i]

  //  Perform butterfly
  __m256d c = _mm256_add_pd(a,b); // c = a + b
  __m256d d = _mm256_sub_pd(a,b); // d = a - b

  _mm256_store_pd((double*)T0, c); // T[c] <- c

  //  Multiply by twiddle factor.
  //  FIXME: Documentation. Tomorrow will have forgotten how this works!
  c = _mm256_mul_pd(d,r);                             // c = (a-b)*omega_r
  d = _mm256_mul_pd(_mm256_shuffle_pd(d,d,0x5), i);   // shuffle -> switch real and imag
  c = _mm256_addsub_pd(c,d);

  _mm256_store_pd((double*)T1, c); // T[c+N/2] <- c
}

// Runs a inverse butterfly on 2 elements simultaniously
// Arguments:
// - twiddles: Vektor of twiddles [r0, i0, r1, i1]
// - T0:       Pointer to 2 complex doubles (first half)
// - T1:       Pointer to 2 complex doubles (second half)
static inline void fft_inverse_butterfly(__m256d twiddle, __m256d* T0, __m256d* T1) {
  __m256d r   = _mm256_unpacklo_pd(twiddle, twiddle);     //   r = [r0,r0,r1,r1]
  __m256d i   = _mm256_unpackhi_pd(twiddle, twiddle);     //   i = [i0,i0,i1,i1]
  i = _mm256_xor_pd(i,_mm256_set1_pd(-0.0));              //   i = -i

  //  Grab elements
  __m256d a = _mm256_load_pd((double*)T0);  // a = [a0r, a0i, a1r, a1i]
  __m256d b = _mm256_load_pd((double*)T1);  // b = [b0r, b0i, b1r, b1i]

  //  Perform butterfly
  __m256d c, d;

  //  Multiply by twiddle factor.
  c = _mm256_mul_pd(b,r);
  d = _mm256_mul_pd(_mm256_shuffle_pd(b,b,0x5),i);
  c = _mm256_addsub_pd(c,d);

  b = _mm256_add_pd(a,c);
  d = _mm256_sub_pd(a,c);

  _mm256_store_pd((double*)T0, b);
  _mm256_store_pd((double*)T1, d);
}

static inline void dft_2p(__m128d* T) {
  __m128d a = T[0];
  __m128d b = T[1];
  T[0] = _mm_add_pd(a,b);
  T[1] = _mm_sub_pd(a,b);
}

void fft_forward(__m128d *T,int k,int tds){
  // (Bit reversed) 2-point DFT
  if(k==1) {
    dft_2p(T);
    return;
  }

  //  Don't thread if it's too small.
  if (k < FFT_THRESHOLD_K) {
    tds = 1;
  }

  size_t length      = 1 << k;
  size_t half_length = length / 2;

  //  Get local twiddle table.
  my_complex_t* local_table = twiddle_table[k];

  for (size_t n = 0; n < half_length; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddles = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    fft_forward_butterfly(twiddles, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }

  if (tds < 2){
    //  No more threads.
    fft_forward(T, k - 1, 1);
    fft_forward(T + half_length, k - 1, 1);
  } else {
    //  Run sub-recursions in parallel.
    int tds0 = tds / 2;
    int tds1 = tds - tds0;
    #pragma omp parallel num_threads(2)
    {
      int tid = omp_get_thread_num();
      if(tid == 0){
          fft_forward(T,k - 1,tds0);
      }
      if(tid != 0 || omp_get_num_threads() < 2){
          fft_forward(T + half_length,k - 1,tds1);
      }
    }
  }
}

// Runs FFT without cached twiddles
// FIXME: Merge?
void fft_forward_uncached(__m128d *T,int k,int tds){
  //  End recursion at 2 points.
  if (k == 1){
    dft_2p(T);
    return;
  }

  //  Don't thread if it's too small.
  if (k < FFT_THRESHOLD_K)
    tds = 1;

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
    __m128d a0 = T[c];
    __m128d b0 = T[c + half_length];

    //  Perform butterfly
    __m128d c0,d0;
    c0 = _mm_add_pd(a0,b0);
    d0 = _mm_sub_pd(a0,b0);

    T[c] = c0;

    //  Multiply by twiddle factor.
    c0 = _mm_mul_pd(d0,r0);
    d0 = _mm_mul_pd(_mm_shuffle_pd(d0,d0,1),i0);
    c0 = _mm_addsub_pd(c0,d0);

    T[c + half_length] = c0;
  }

  if (tds < 2){
    //  No more threads.
    if(k-1 < twiddle_table_size) {
      fft_forward(T, k - 1, 1);
      fft_forward(T + half_length, k - 1, 1);
    } else {
      fft_forward_uncached(T, k - 1, 1);
      fft_forward_uncached(T + half_length, k - 1, 1);
    }
  }else{
    //  Run sub-recursions in parallel.
    int tds0 = tds / 2;
    int tds1 = tds - tds0;
    #pragma omp parallel num_threads(2)
    {
      int tid = omp_get_thread_num();
      if(tid == 0){
        if(k-1 < twiddle_table_size) {
          fft_forward(T,k - 1,tds0);
        } else {
          fft_forward_uncached(T,k - 1,tds0);
        }
      }
      if(tid != 0 || omp_get_num_threads() < 2){
        if(k-1 < twiddle_table_size) {
          fft_forward(T + half_length,k - 1,tds1);
        } else {
          fft_forward_uncached(T + half_length,k - 1,tds1);
        }
      }
    }
  }
}

void fft_inverse(__m128d *T,int k,int tds){
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

  //  Don't thread if it's too small.
  if (k < FFT_THRESHOLD_K)
      tds = 1;

  size_t length = 1 << k;
  size_t half_length = length / 2;

  if (tds < 2){
      //  No more threads.
      fft_inverse(T,k - 1,1);
      fft_inverse(T + half_length,k - 1,1);
  }else{
      //  Run sub-recursions in parallel.
      int tds0 = tds / 2;
      int tds1 = tds - tds0;
      #pragma omp parallel num_threads(2)
      {
          int tid = omp_get_thread_num();
          if (tid == 0){
              fft_inverse(T,k - 1,tds0);
          }
          if (tid != 0 || omp_get_num_threads() < 2){
              fft_inverse(T + half_length,k - 1,tds1);
          }
      }
  }

  //  Get local twiddle table.
  my_complex_t* local_table = twiddle_table[k];

  //  Perform FFT reduction into two halves.
  for (size_t n = 0; n < half_length; n+=2){
      //  Grab Twiddle Factors
      __m256d twiddle = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
      fft_inverse_butterfly(twiddle, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }
}

void fft_inverse_uncached(__m128d *T,int k,int tds){
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

  //  Don't thread if it's too small.
  if (k < FFT_THRESHOLD_K)
      tds = 1;

  size_t length = 1 << k;
  double omega = 2 * M_PI / length;
  size_t half_length = length / 2;

  if (tds < 2){
      //  No more threads.
      if(k - 1 < twiddle_table_size) {
        fft_inverse(T,k - 1,1);
        fft_inverse(T + half_length,k - 1,1);
      } else {
        fft_inverse_uncached(T,k - 1,1);
        fft_inverse_uncached(T + half_length,k - 1,1);
      }
  }else{
      //  Run sub-recursions in parallel.
      int tds0 = tds / 2;
      int tds1 = tds - tds0;
      #pragma omp parallel num_threads(2)
      {
          int tid = omp_get_thread_num();
          if (tid == 0){
            if(k - 1 < twiddle_table_size) {
              fft_inverse(T,k - 1,tds0);
            } else {
              fft_inverse_uncached(T,k - 1,tds0);
            }
          }
          if (tid != 0 || omp_get_num_threads() < 2){
            if(k - 1 < twiddle_table_size) {
              fft_inverse(T + half_length,k - 1,tds1);
            } else {
              fft_inverse_uncached(T + half_length,k - 1,tds1);
            }
          }
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
      __m128d a0 = T[c];
      __m128d b0 = T[c + half_length];

      //  Perform butterfly
      __m128d c0,d0;

      //  Multiply by twiddle factor.
      c0 = _mm_mul_pd(b0,r0);
      d0 = _mm_mul_pd(_mm_shuffle_pd(b0,b0,1),i0);
      c0 = _mm_addsub_pd(c0,d0);

      b0 = _mm_add_pd(a0,c0);
      d0 = _mm_sub_pd(a0,c0);

      T[c] = b0;
      T[c + half_length] = d0;
  }
}

