#include <stdint.h>
#include <stdio.h>
#include "fft.h"

// Puts 12 bits per complex point
//
// |                word 1                |               word 2
// | ******** **** **** ******** ******** | **** **** ******** ******** **** .... --->
// |              |             |               |             |             |
// |      T[0]    |     T[1]    |     T[2]      |     T[3]    |     T[4]    |
//
//      word 2     |                word 3                |
// <--- .... ***** | ******** ******** **** **** ******** |
//           |                |             |              |
//           |      T[5]      |    T[6]     |     T[7]     |
static inline size_t int_to_fft12(complex double *V, const uint32_t *A, size_t AL) {
  __m128d* T = (__m128d*)V;
  complex double *origV = V;

  for (size_t c = 0; c < AL/3+2; c++){
    if(3*c >= AL) {
      break;
    }

    uint32_t word1 = A[3*c];

    *T++ = _mm_set_sd(word1 & 0xfff);
    word1 <<= 12;
    *T++ = _mm_set_sd(word1 & 0xfff);
    word1 <<= 12;

    if(3*c+1 >= AL) {
      *T++ = _mm_set_sd(word1 & 0xff);
      break;
    }

    uint32_t word2 = A[3*c+1];

    *T++ = _mm_set_sd((word1 & 0xff) | (word2 & 0xf) >> 8);
    word2 <<= 4;

    *T++ = _mm_set_sd(word2 & 0xfff);
    word2 <<= 12;

    *T++ = _mm_set_sd(word2 & 0xfff);
    word2 <<= 12;

    if(3*c+2 >= AL) {
      *T++ = _mm_set_sd(word2 & 0xf);
      break;
    }

    uint32_t word3 = A[3*c+2];

    *T++ = _mm_set_sd((word2 & 0xf) | (word3 & 0xff) >> 8);
    word3 <<= 8;

    *T++ = _mm_set_sd(word2 & 0xfff);
    word3 <<= 12;

    *T++ = _mm_set_sd(word2 & 0xfff);
    word3 <<= 12;
  }

  V = (complex double*)T;
  return V-origV;
}

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
    word1 += (carry % 256) * 256 * 256;     //  Get 8 bits.
    carry /= 256;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 += (carry % 256) * 256 * 256 * 256;  //  Get 8 bits.
    carry /= 256;

    A[c] = word1;
  }
}

static inline void fft_to_int12(const complex double *T, uint32_t *A, size_t AL, double scale) {
  uint64_t carry = 0;

  for (size_t c = 0; c < AL/3+2; c++){
    double   f_point;
    uint64_t i_point;
    uint32_t word1, word2, word3;

    if(3*c >= AL) {
      break;
    }

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 = carry & 0xfff;                  //  Get 12 bits.
    carry <<= 12;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 |= (carry & 0xfff) << 12;         //  Get 12 bits.
    carry <<= 12;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word1 |= (carry & 0xff) << 24;          //  Get 8+4 bits.
    carry <<= 8;
    A[3*c] = word1;
    if(3*c+1 >= AL) {
      break;
    }

    word2 = carry & 0xf;                    // FIXME: DETECT END!
    carry <<= 4;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word2 |= (carry & 0xfff) << 4;          //  Get 12 bits.
    carry <<= 12;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word2 |= (carry & 0xfff) << 16;         //  Get 12 bits.
    carry <<= 12;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word2 |= (carry & 0xf) << 28;           //  Get 4+8 bits.
    carry <<= 4;
    A[3*c+1] = word2;
    if(3*c+2 >= AL) {
      break;
    }

    word3 = carry & 0xff;                   // FIXME: DETECT END!
    carry <<= 8;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word3 |= (carry & 0xfff) << 8;          //  Get 12 bits.
    carry <<= 12;

    f_point = ((double*)T++)[0] * scale;    //  Load and scale
    i_point = (uint64_t)(f_point + 0.5);    //  Round
    carry += i_point;                       //  Add to carry
    word3 |= (carry & 0xfff) << 20;         //  Get 12 bits.
    A[3*c+2] = word3;
  }
}

// Converts an array of words to an array of complex numbers. Put a given
// number of digits per point.
size_t int_to_fft(complex double *V, int k, const uint32_t *A, size_t AL, int bits_per_point) {
  size_t points_written = 0;

  switch(bits_per_point) {
    case 8:  points_written = int_to_fft8(V, A, AL);  break;
    case 12: points_written = int_to_fft12(V, A, AL); break;
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
    case 8:  fft_to_int8(T, A, AL, scale);  break;
    case 12: fft_to_int12(T, A, AL, scale); break;
    default:
      fprintf(stderr, "Not implemented\n");
      abort();
    break;
  }
}
