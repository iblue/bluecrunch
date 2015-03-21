#include <stdint.h>
#include <x86intrin.h>
#include "fft.h"

// Puts 2 digits per complex point
static inline size_t int_to_fft2(complex double *V, const uint32_t *A, size_t AL) {
  __m128d* T = (__m128d*)V;
  complex double *origV = V;

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

  V = (complex double*)T;
  return V-origV;
}

// Puts 3 digits per complex point
static inline size_t int_to_fft3(complex double *V, const uint32_t *A, size_t AL) {
  __m128d* T = (__m128d*)V;
  complex double *origV = V;


  for (size_t c = 0; c < AL; c++) {
    uint32_t word = A[c];

    *T++ = _mm_set_sd(word % 1000);
    word /= 1000;
    *T++ = _mm_set_sd(word % 1000);
    word /= 1000;
    *T++ = _mm_set_sd(word);
  }

  V = (complex double*)T;
  return V-origV;
}

// Puts 4 digits per complex point
static inline size_t int_to_fft4(complex double *V, const uint32_t *A, size_t AL) {
  __m128d* T = (__m128d*)V;
  complex double *origV = V;

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

  V = (complex double*)T;
  return V-origV;
}

static inline void fft_to_int2(const complex double *T, uint32_t *A, size_t AL, double scale) {
  uint64_t carry = 0;

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

static inline void fft_to_int3(const complex double *T, uint32_t *A, size_t AL, double scale) {
  uint64_t carry = 0;

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

static inline void fft_to_int4(const complex double *T, uint32_t *A, size_t AL, double scale) {
  uint64_t carry = 0;

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
  }
}

// Converts an array of words to an array of complex numbers. Put a given
// number of digits per point.
size_t int_to_fft(complex double *V, int k, const uint32_t *A, size_t AL, int digits_per_point) {
  size_t points_written = 0;

  switch(digits_per_point) {
    case 2: points_written = int_to_fft2(V, A, AL); break;
    case 3: points_written = int_to_fft3(V, A, AL); break;
    case 4: points_written = int_to_fft4(V, A, AL); break;
  }

  //  Pad the rest with zeros.
  __m128d* T = (__m128d*)V;
  size_t fft_length = 1 << k;

  for(size_t i = points_written; i < fft_length; i++) {
    T[i] = _mm_setzero_pd();
  }

  return points_written;
}

//  Convert FFT array back to word array. Perform rounding and carryout.
void fft_to_int(const complex double *T, int k, uint32_t *A, size_t AL, int digits_per_point) {
  //  Compute Scaling Factor
  size_t fft_length = 1 << k;
  double scale = 1. / fft_length;

  //  Round and carry out.
  switch(digits_per_point) {
    case 2: fft_to_int2(T, A, AL, scale); break;
    case 3: fft_to_int3(T, A, AL, scale); break;
    case 4: fft_to_int4(T, A, AL, scale); break;
  }
}
