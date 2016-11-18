#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> // abort
#include "fft.h"

#include "experiments/codegen/fft_to_int7.c"
#include "experiments/codegen/fft_to_int8.c"
#include "experiments/codegen/fft_to_int9.c"
#include "experiments/codegen/fft_to_int10.c"
#include "experiments/codegen/fft_to_int11.c"
#include "experiments/codegen/fft_to_int12.c"
#include "experiments/codegen/fft_to_int13.c"
#include "experiments/codegen/fft_to_int14.c"
#include "experiments/codegen/fft_to_int15.c"
#include "experiments/codegen/fft_to_int16.c"
#include "experiments/codegen/fft_to_int17.c"
#include "experiments/codegen/fft_to_int18.c"
#include "experiments/codegen/fft_to_int19.c"
#include "experiments/codegen/fft_to_int20.c"
#include "experiments/codegen/fft_to_int21.c"
#include "experiments/codegen/fft_to_int22.c"

#include "experiments/codegen/int_to_fft7.c"
#include "experiments/codegen/int_to_fft8.c"
#include "experiments/codegen/int_to_fft9.c"
#include "experiments/codegen/int_to_fft10.c"
#include "experiments/codegen/int_to_fft11.c"
#include "experiments/codegen/int_to_fft12.c"
#include "experiments/codegen/int_to_fft13.c"
#include "experiments/codegen/int_to_fft14.c"
#include "experiments/codegen/int_to_fft15.c"
#include "experiments/codegen/int_to_fft16.c"
#include "experiments/codegen/int_to_fft17.c"
#include "experiments/codegen/int_to_fft18.c"
#include "experiments/codegen/int_to_fft19.c"
#include "experiments/codegen/int_to_fft20.c"
#include "experiments/codegen/int_to_fft21.c"
#include "experiments/codegen/int_to_fft22.c"


// Converts an array of words to an array of complex numbers. Put a given
// number of digits per point.
size_t int_to_fft(complex double *V, size_t length, const uint32_t *A, size_t AL, int bits_per_point) {
  size_t points_written = 0;

  switch(bits_per_point) {
    case 7:  points_written = int_to_fft7(V, A, AL);  break;
    case 8:  points_written = int_to_fft8(V, A, AL);  break;
    case 9:  points_written = int_to_fft9(V, A, AL);  break;
    case 10: points_written = int_to_fft10(V, A, AL); break;
    case 11: points_written = int_to_fft11(V, A, AL); break;
    case 12: points_written = int_to_fft12(V, A, AL); break;
    case 13: points_written = int_to_fft13(V, A, AL); break;
    case 14: points_written = int_to_fft14(V, A, AL); break;
    case 15: points_written = int_to_fft15(V, A, AL); break;
    case 16: points_written = int_to_fft16(V, A, AL); break;
    case 17: points_written = int_to_fft17(V, A, AL); break;
    case 18: points_written = int_to_fft18(V, A, AL); break;
    case 19: points_written = int_to_fft19(V, A, AL); break;
    case 20: points_written = int_to_fft20(V, A, AL); break;
    case 21: points_written = int_to_fft21(V, A, AL); break;
    case 22: points_written = int_to_fft22(V, A, AL); break;
    default:
      fprintf(stderr, "Not implemented\n");
      abort();
    break;
  }

  //  Pad the rest with zeros.
  __m128d* T = (__m128d*)V;

  for(size_t i = points_written; i < length; i++) {
    T[i] = _mm_setzero_pd();
  }

  return points_written;
}

//  Convert FFT array back to word array. Perform rounding and carryout.
void fft_to_int(const complex double *T, size_t length, uint32_t *A, size_t AL, int bits_per_point) {
  //  Compute Scaling Factor
  double scale = 1. / length;

  //  Round and carry out.
  switch(bits_per_point) {
    case 7:  fft_to_int7(T, A, AL, scale);  break;
    case 8:  fft_to_int8(T, A, AL, scale);  break;
    case 9:  fft_to_int9(T, A, AL, scale);  break;
    case 10: fft_to_int10(T, A, AL, scale); break;
    case 11: fft_to_int11(T, A, AL, scale); break;
    case 12: fft_to_int12(T, A, AL, scale); break;
    case 13: fft_to_int13(T, A, AL, scale); break;
    case 14: fft_to_int14(T, A, AL, scale); break;
    case 15: fft_to_int15(T, A, AL, scale); break;
    case 16: fft_to_int16(T, A, AL, scale); break;
    case 17: fft_to_int17(T, A, AL, scale); break;
    case 18: fft_to_int18(T, A, AL, scale); break;
    case 19: fft_to_int19(T, A, AL, scale); break;
    case 20: fft_to_int20(T, A, AL, scale); break;
    case 21: fft_to_int21(T, A, AL, scale); break;
    case 22: fft_to_int22(T, A, AL, scale); break;
    default:
      fprintf(stderr, "Not implemented\n");
      abort();
    break;
  }
}
