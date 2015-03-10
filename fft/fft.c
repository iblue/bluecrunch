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

my_complex_t* twiddle_table[32];
int twiddle_table_size = 0;

void fft_ensure_table(int k) {
  //  Makes sure the twiddle factor table is large enough to handle an FFT of
  //  size 2^k.

  //  Do one level at a time
  if (k - 1 > 0) {
    fft_ensure_table(k - 1);
  }

  size_t length = 1 << k;
  double omega = 2 * M_PI / length;
  length /= 2;

  //  Build the sub-table.
  my_complex_t *sub_table = (my_complex_t*)malloc(length*sizeof(my_complex_t));

  for (size_t c = 0; c < length; c++){
      //  Generate Twiddle Factor
      double angle = omega * c;
      my_complex_t twiddle_factor;
      twiddle_factor.r = cos(angle);
      twiddle_factor.i = sin(angle);
      sub_table[c] = twiddle_factor;
  }

  //  Push into main table.
  twiddle_table[k] = sub_table;
  twiddle_table_size = k+1;
}

void fft4(__m128d *T) {
  // Bit reversed 4-point DFT
  // T[0] = a +   b + c +   d (= x_0)
  // T[1] = a -   b + c -   d (= x_2)
  // T[2] = a + i*b + c + i*d (= x_1)
  // T[3] = a - i*b - c + i*d (= x_3)
  //
  // Let [0] bet the real, [1] be the complex part:
  // A[0] = a[0] + b[0] + c[0] + d[0]; A[1] = a[1] + b[1] + c[1] + d[1];
  // B[0] = a[0] - b[0] + c[0] - d[0]; B[1] = a[1] - b[1] + c[1] - d[1];
  // C[0] = a[0] - b[1] - c[0] + d[1]; C[1] = a[1] + b[0] - c[1] - d[0];
  // D[0] = a[0] + b[1] - c[0] - d[1]; D[1] = a[1] - b[0] - c[1] + d[0];

  // Input
  double a[2] = {((double*)&T[0])[0], ((double*)&T[0])[1]};
  double b[2] = {((double*)&T[1])[0], ((double*)&T[1])[1]};
  double c[2] = {((double*)&T[2])[0], ((double*)&T[2])[1]};
  double d[2] = {((double*)&T[3])[0], ((double*)&T[3])[1]};

  // Calc
  double A[2], B[2], C[2], D[2];

  A[0] = a[0] + b[0] + c[0] + d[0]; A[1] = a[1] + b[1] + c[1] + d[1];
  B[0] = a[0] - b[0] + c[0] - d[0]; B[1] = a[1] - b[1] + c[1] - d[1];
  C[0] = a[0] - b[1] - c[0] + d[1]; C[1] = a[1] + b[0] - c[1] - d[0];
  D[0] = a[0] + b[1] - c[0] - d[1]; D[1] = a[1] - b[0] - c[1] + d[0];

  // Output
  T[0] = (__m128d){A[0], A[1]};
  T[1] = (__m128d){B[0], B[1]};
  T[2] = (__m128d){C[0], C[1]};
  T[3] = (__m128d){D[0], D[1]};
}

void fft_forward(__m128d *T,int k,int tds){
  if (k == 2) {
    fft4(T);

    return;
  }

  // (Bit reversed) 2-point DFT
  if (k == 1){
    __m128d a = T[0];
    __m128d b = T[1];
    T[0] = _mm_add_pd(a,b);
    T[1] = _mm_sub_pd(a,b);
    return;
  }

  //  Don't thread if it's too small.
  if (k < FFT_THRESHOLD_K)
    tds = 1;

  size_t length = 1 << k;
  size_t half_length = length / 2;

  //  Get local twiddle table.
  my_complex_t* local_table = twiddle_table[k];

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
  for (size_t n = 0; n < half_length; n+=2){
    //  Grab Twiddle Factors
    __m256d tmp = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    __m256d r   = _mm256_unpacklo_pd(tmp, tmp);             //   r = [r0,r0,r1,r1]
    __m256d i   = _mm256_unpackhi_pd(tmp, tmp);             //   i = [i0,i0,i1,i1]

    //  Grab elements
    __m256d a = _mm256_load_pd((double*)&T[n]);              // a = [a0r, a0i, a1r, a1i]
    __m256d b = _mm256_load_pd((double*)&T[n+half_length]);  // b = [b0r, b0i, b1r, b1i]

    //  Check: Grab Twiddle Factor
    __m128d r0 = _mm_loaddup_pd(&local_table[n].r);
    __m128d i0 = _mm_loaddup_pd(&local_table[n].i);
    __m128d r1 = _mm_loaddup_pd(&local_table[n+1].r);
    __m128d i1 = _mm_loaddup_pd(&local_table[n+1].i);

    //  Check: Grab elements
    __m128d a0 = T[n];
    __m128d b0 = T[n + half_length];
    __m128d a1 = T[n+1];
    __m128d b1 = T[n+1 + half_length];

    //  Perform butterfly
    __m256d c = _mm256_add_pd(a,b); // c = a + b
    __m256d d = _mm256_sub_pd(a,b); // d = a - b

    //  Check: Perform butterfly
    __m128d c0,d0;
    c0 = _mm_add_pd(a0,b0);
    d0 = _mm_sub_pd(a0,b0);
    __m128d c1,d1;
    c1 = _mm_add_pd(a1,b1);
    d1 = _mm_sub_pd(a1,b1);

    _mm256_store_pd((double*)(T+n), c); // T[c] <- c
    T[n] = c0;
    T[n+1] = c1;

    //  Multiply by twiddle factor.
    c = _mm256_mul_pd(d,r);
    d = _mm256_mul_pd(_mm256_shuffle_pd(d,d,1),i);
    c = _mm256_addsub_pd(c,d);

    //  Check: Multiply by twiddle factor.
    c0 = _mm_mul_pd(d0,r0);
    d0 = _mm_mul_pd(_mm_shuffle_pd(d0,d0,1),i0);
    c0 = _mm_addsub_pd(c0,d0);
    c1 = _mm_mul_pd(d1,r1);
    d1 = _mm_mul_pd(_mm_shuffle_pd(d1,d1,1),i1);
    c1 = _mm_addsub_pd(c1,d1);

    _mm256_store_pd((double*)(T+half_length+n), c); // T[c+N/2] <- c
    T[n+half_length] = c0;
    T[n+half_length+1] = c1;
  }

  if (tds < 2){
    //  No more threads.
    fft_forward(T, k - 1, 1);
    fft_forward(T + half_length, k - 1, 1);
  }else{
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

void fft_forward_uncached(__m128d *T,int k,int tds){
  //  Fast Fourier Transform
  //  This function performs a forward FFT of length 2^k.

  //  This is a Decimation-in-Frequency (DIF) FFT.
  //  The frequency domain output is in bit-reversed order.

  //Parameters:
  //  -   T           -   Pointer to array.
  //  -   k           -   2^k is the size of the transform

  //  End recursion at 2 points.
  if (k == 1){
    __m128d a = T[0];
    __m128d b = T[1];
    T[0] = _mm_add_pd(a,b);
    T[1] = _mm_sub_pd(a,b);
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
      __m128d a = T[0];
      __m128d b = T[1];
      T[0] = _mm_add_pd(a,b);
      T[1] = _mm_sub_pd(a,b);
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
  for (size_t c = 0; c < half_length; c++){
      //  Grab Twiddle Factor
      __m128d r0 = _mm_loaddup_pd(&local_table[c].r);
      __m128d i0 = _mm_loaddup_pd(&local_table[c].i);
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
      __m128d a = T[0];
      __m128d b = T[1];
      T[0] = _mm_add_pd(a,b);
      T[1] = _mm_sub_pd(a,b);
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

void fft_pointwise(__m128d *T,__m128d *A,int k){
  //  Performs pointwise multiplications of two FFT arrays.

  //Parameters:
  //  -   T           -   Pointer to array.
  //  -   k           -   2^k is the size of the transform

  size_t length = 1 << k;
  for (size_t c = 0; c < length; c++){
      __m128d a0 = T[c];
      __m128d b0 = A[c];
      __m128d c0,d0;
      c0 = _mm_mul_pd(a0,_mm_unpacklo_pd(b0,b0));
      d0 = _mm_mul_pd(_mm_shuffle_pd(a0,a0,1),_mm_unpackhi_pd(b0,b0));
      T[c] = _mm_addsub_pd(c0,d0);
  }
}
void int_to_fft(__m128d *T,int k,const uint32_t *A,size_t AL, int digits_per_point){
  //  Convert word array into FFT array. Put 2 decimal digits per complex point.

  //Parameters:
  //  -   T   -   FFT array
  //  -   k   -   2^k is the size of the transform
  //  -   A   -   word array
  //  -   AL  -   length of word array

  size_t fft_length = 1 << k;
  __m128d *Tstop = T + fft_length;

  //  Since there are 9 digits per word and we want to put 2 digits per
  //  point, the length of the transform must be at least 4.5 times the word
  //  length of the input.
  if (fft_length < (9/digits_per_point)*AL) {
    abort();
  }

  //  Convert
  if(digits_per_point == 2) {
    for (size_t c = 0; c < AL/2+1; c++){
      uint32_t word1 = A[2*c];

      if(2*c >= AL) {
        break;
      }

      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;

      if(2*c+1 >= AL) {
        *T++ = _mm_set_sd(word1);
        break;
      }

      uint32_t word2 = A[2*c+1];
      uint32_t tmp = word2 % 10 * 10 + word1;
      *T++ = _mm_set_sd(tmp);

      word2 /= 10;

      *T++ = _mm_set_sd(word2 % 100);
      word2 /= 100;
      *T++ = _mm_set_sd(word2 % 100);
      word2 /= 100;
      *T++ = _mm_set_sd(word2 % 100);
      word2 /= 100;
      *T++ = _mm_set_sd(word2 % 100);
      word2 /= 100;
    }
  }

  if(digits_per_point == 3) {
    for (size_t c = 0; c < AL; c++) {
      uint32_t word = A[c];

      *T++ = _mm_set_sd(word % 1000);
      word /= 1000;
      *T++ = _mm_set_sd(word % 1000);
      word /= 1000;
      *T++ = _mm_set_sd(word);
    }
  }

  if(digits_per_point == 4) {
    for (size_t c = 0; c < AL/4+1; c++){
      uint32_t tmp1, tmp2, tmp3;

      if(4*c >= AL) {
        break;
      }

      uint32_t word1 = A[4*c];
      *T++ = _mm_set_sd(word1 % 10000);   // W1: 4
      word1 /= 10000;
      *T++ = _mm_set_sd(word1 % 10000);   // W1: 4
      word1 /= 10000;

      if(4*c+1 >= AL) {
        *T++ = _mm_set_sd(word1);         // branch out: W1: 1
        break;
      }

      uint32_t word2 = A[4*c+1];
      tmp1 = word2 % 1000 * 10 + word1;   // W1: 1, W2: 3
      *T++ = _mm_set_sd(tmp1);
      word2 /= 1000;
      *T++ = _mm_set_sd(word2 % 10000);   // W2: 4
      word2 /= 10000;


      if(4*c+2 >= AL) {
        *T++ = _mm_set_sd(word2);         // branch out: W2: 2
        break;
      }

      uint32_t word3 = A[4*c+2];
      tmp2 = word3 % 100 * 100 + word2;
      *T++ = _mm_set_sd(tmp2);            // W2: 2, W3: 2
      word3 /= 100;
      *T++ = _mm_set_sd(word3 % 10000);   // W3: 4
      word3 /= 10000;

      if(4*c+3 >= AL) {
        *T++ = _mm_set_sd(word3);         // branch out: W3: 3
        break;
      }

      uint32_t word4 = A[4*c+3];
      tmp3 = word4 % 10 * 1000 + word3;
      *T++ = _mm_set_sd(tmp3);            // W3: 3, W4: 1
      word4 /= 10;
      *T++ = _mm_set_sd(word4 % 10000);   // W4: 4
      word4 /= 10000;
      *T++ = _mm_set_sd(word4 % 10000);   // W4: 4
      word4 /= 10000;
    }
  }

  //  Pad the rest with zeros.
  while (T < Tstop)
    *T++ = _mm_setzero_pd();
}

void fft_to_int(__m128d *T,int k,uint32_t *A,size_t AL, int digits_per_point){
  //  Convert FFT array back to word array. Perform rounding and carryout.

  //Parameters:
  //  -   T   -   FFT array
  //  -   A   -   word array
  //  -   AL  -   length of word array

  //  Compute Scaling Factor
  size_t fft_length = 1 << k;
  double scale = 1. / fft_length;

  //  Since there are 9 digits per word and we want to put 3 digits per
  //  point, the length of the transform must be at least 3 times the word
  //  length of the input.
  if (fft_length < (9/digits_per_point)*AL) {
    abort();
  }

  //  Round and carry out.
  uint64_t carry = 0;
  if(digits_per_point == 2) {
    for (size_t c = 0; c < AL/2+1; c++){
      double   f_point;
      uint64_t i_point;
      uint32_t word1;
      uint32_t word2;

      if(2*c >= AL) {
        break;
      }

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 = carry % 100;                    //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 100) * 100;           //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 100) * 10000;         //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 100) * 1000000;       //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 10) * 100000000;      //  Get 1 digit.
      carry /= 10;

      A[2*c] = word1;

      if(2*c+1 >= AL) {
        break;
      }

      word2  = carry % 10;                    // Get 1 digit
      carry /= 10;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 100) * 10;            //  Get 2 digits
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 100) * 1000;         //  Get 2 digits
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 100) * 100000;        //  Get 2 digits
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 100) * 10000000;      //  Get 2 digits
      carry /= 100;

      A[2*c+1] = word2;
    }
  }

  if(digits_per_point == 3) {
    for (size_t c = 0; c < AL; c++){
        double   f_point;
        uint64_t i_point;
        uint32_t word;

        f_point = ((double*)T++)[0] * scale;    //  Load and scale
        i_point = (uint64_t)(f_point + 0.5);    //  Round
        carry += i_point;                       //  Add to carry
        word = carry % 1000;                    //  Get 3 digits.
        carry /= 1000;

        f_point = ((double*)T++)[0] * scale;    //  Load and scale
        i_point = (uint64_t)(f_point + 0.5);    //  Round
        carry += i_point;                       //  Add to carry
        word += (carry % 1000) * 1000;          //  Get 3 digits.
        carry /= 1000;

        f_point = ((double*)T++)[0] * scale;    //  Load and scale
        i_point = (uint64_t)(f_point + 0.5);    //  Round
        carry += i_point;                       //  Add to carry
        word += (carry % 1000) * 1000000;       //  Get 3 digits.
        carry /= 1000;

        A[c] = word;
    }
  }

  if(digits_per_point == 4) {
    for (size_t c = 0; c < AL/4+1; c++){
      double   f_point;
      uint64_t i_point;
      uint32_t word1, word2, word3, word4;

      if(4*c >= AL) {
        break;
      }

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 = carry % 10000;                  //  Get 4 digits
      carry /= 10000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 10000) * 10000;       //  Get 4 digits
      carry /= 10000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word1 += (carry % 10) * 100000000;      //  Get 1 digit.
      carry /= 10;

      A[4*c]   = word1;

      if(word1 >= 1000000000) {
        abort();
      }

      if(4*c+1 >= AL) {
        break;
      }

      word2  = carry % 1000;                  //  Get 3 digits.
      carry /= 1000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 10000) * 1000;        //  Get 4 digits
      carry /= 10000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word2 += (carry % 100) * 10000000;      //  Get 2 digits
      carry /= 100;

      A[4*c+1] = word2;

      if(word2 >= 1000000000) {
        abort();
      }

      if(4*c+2 >= AL) {
        break;
      }

      word3  = carry % 100;                   //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word3 += (carry % 10000) * 100;         //  Get 4 digits
      carry /= 10000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word3 += (carry % 1000) * 1000000;      //  Get 3 digits
      carry /= 1000;

      A[4*c+2] = word3;

      if(word3 >= 1000000000) {
        abort();
      }

      if(4*c+3 >= AL) {
        break;
      }

      word4  = carry % 10;                    //  Get 1 digit
      carry /= 10;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word4 += (carry % 10000) * 10;          //  Get 4 digits
      carry /= 10000;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word4 += (carry % 10000) * 100000;      //  Get 4 digits
      carry /= 10000;

      A[4*c+3] = word4;

      if(word4 >= 1000000000) {
        abort();
      }
    }
  }
}
