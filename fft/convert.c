#include <stdint.h>
#include <stdio.h>
#include "fft.h"

// Puts 8 bits per complex point
static inline size_t int_to_fft8(complex double *V, const uint32_t *A, size_t AL) {
  __m128d* T = (__m128d*)V;
  complex double *origV = V;

  for (size_t c = 0; c < AL; c++){
    uint32_t word = A[c];

    *T++ = _mm_set_sd(word % 256);
    word /= 256;
    *T++ = _mm_set_sd(word % 256);
    word /= 256;
    *T++ = _mm_set_sd(word % 256);
    word /= 256;
    *T++ = _mm_set_sd(word % 256);
    word /= 256;
  }

  V = (complex double*)T;
  return V-origV;
}

static inline void fft_to_int8(const complex double *T, uint32_t *A, size_t AL, double scale) {
  uint64_t carry = 0;

  for (size_t c = 0; c < AL; c++){
    double   f_point;
    uint64_t i_point;
    uint32_t word1;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 = carry % 256;                    //  Get 8 bits.
    carry /= 256;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 += (carry % 256) * 256;           //  Get 8 bits.
    carry /= 256;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 += (carry % 100) * 256 * 256;     //  Get 8 bits.
    carry /= 256;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 += (carry % 100) * 256 * 256 * 256;  //  Get 8 bits.
    carry /= 256;

    A[c] = word1;
  }
}

// Converts an array of words to an array of complex numbers. Put a given
// number of digits per point.
size_t int_to_fft(complex double *V, int k, const uint32_t *A, size_t AL, int bits_per_point) {
  size_t points_written = 0;

  switch(bits_per_point) {
    case 8: points_written = int_to_fft8(V, A, AL); break;
    default:
      fprintf(stderr, "Not implemented\n");
      abort();
    break;
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
void fft_to_int(const complex double *T, int k, uint32_t *A, size_t AL, int bits_per_point) {
  //  Compute Scaling Factor
  size_t fft_length = 1 << k;
  double scale = 1. / fft_length;

  //  Round and carry out.
  switch(bits_per_point) {
    case 8: fft_to_int8(T, A, AL, scale); break;
    default:
      fprintf(stderr, "Not implemented\n");
      abort();
    break;
  }
}
