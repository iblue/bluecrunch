#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <malloc.h>
#include <pmmintrin.h>
#include <omp.h>
#include "fft.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

namespace Bluecrunch {
  std::vector<std::vector<my_complex>> twiddle_table;

  void fft_ensure_table(int k){
    //  Makes sure the twiddle factor table is large enough to handle an FFT of
    //  size 2^k.

    int current_k = (int)twiddle_table.size() - 1;
    if (current_k >= k)
        return;

    //  Do one level at a time
    if (k - 1 > current_k){
        fft_ensure_table(k - 1);
    }

    size_t length = 1 << k;
    double omega = 2 * M_PI / length;
    length /= 2;

    //  Build the sub-table.
    std::vector<my_complex> sub_table;
    for (size_t c = 0; c < length; c++){
        //  Generate Twiddle Factor
        double angle = omega * c;
        auto twiddle_factor = my_complex(cos(angle),sin(angle));
        sub_table.push_back(twiddle_factor);
    }

    //  Push into main table.
    twiddle_table.push_back(std::move(sub_table));
  }

  void fft_forward(__m128d *T,int k,int tds){
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
    size_t half_length = length / 2;

    //  Get local twiddle table.
    std::vector<my_complex> &local_table = twiddle_table[k];

    //  Perform FFT reduction into two halves.
    for (size_t c = 0; c < half_length; c++){
      //  Grab Twiddle Factor
      __m128d r0 = _mm_loaddup_pd(&local_table[c].r);
      __m128d i0 = _mm_loaddup_pd(&local_table[c].i);

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
      fft_forward(T,k - 1);
      fft_forward(T + half_length,k - 1);
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
        fft_inverse(T,k - 1);
        fft_inverse(T + half_length,k - 1);
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
    std::vector<my_complex> &local_table = twiddle_table[k];

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
void int_to_fft(__m128d *T,int k,const uint32_t *A,size_t AL){
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
    if (fft_length < 5*AL)
        throw "FFT length is too small.";

    //  Convert
    for (size_t c = 0; c < AL/2; c++){
      uint32_t word1 = A[2*c];
      uint32_t word2 = A[2*c+1];

      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;
      *T++ = _mm_set_sd(word1 % 100);
      word1 /= 100;


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

    if(AL%2 == 1) {
      uint32_t word = A[AL-1];

      *T++ = _mm_set_sd(word % 100);
      word /= 100;
      *T++ = _mm_set_sd(word % 100);
      word /= 100;
      *T++ = _mm_set_sd(word % 100);
      word /= 100;
      *T++ = _mm_set_sd(word % 100);
      word /= 100;
      *T++ = _mm_set_sd(word);
    }

    //  Pad the rest with zeros.
    while (T < Tstop)
      *T++ = _mm_setzero_pd();
}

void fft_to_int(__m128d *T,int k,uint32_t *A,size_t AL){
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
    if (fft_length < 5*AL)
        throw "FFT length is too small.";

    //  Round and carry out.
    uint64_t carry = 0;
    for (size_t c = 0; c < AL/2; c++){
        double   f_point;
        uint64_t i_point;
        uint32_t word1;
        uint32_t word2;

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

    if(AL%2 == 1) {
      double   f_point;
      uint64_t i_point;
      uint32_t word;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word  = carry % 100;                    //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word  += (carry % 100) * 100;           //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word  += (carry % 100) * 10000;         //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word  += (carry % 100) * 1000000;       //  Get 2 digits.
      carry /= 100;

      f_point = ((double*)T++)[0] * scale;    //  Load and scale
      i_point = (uint64_t)(f_point + 0.5);    //  Round
      carry += i_point;                       //  Add to carry
      word  += (carry % 10) * 100000000;      //  Get 1 digit.
      carry /= 10;

      // FIXME: carry should now == 0

      A[AL-1] = word;
    }
  }
}
