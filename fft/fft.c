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

// When do we start to decompose FFTs for parallization
#define FFT_THRESHOLD_LENGTH 60000

// Find optimal transform length
size_t fft_length(size_t source_length) {
  size_t length = 1;
  while(1) {
    if(2*length >= source_length) {
      return 2*length;
    }
    // 3*2^k transform
    if(3*length >= source_length) {
      return 3*length;
    }
    length *= 2;
  }
}

void _fft_forward(complex double *T, size_t length) {
  // (Bit reversed) 2-point DFT
  if(length == 2) {
    dft_2p(T);
    return;
  }

  if(length == 3) {
    dft_3p(T);
    return;
  }

  size_t half_length = length / 2;

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[table_select(length)];

  if(half_length > 32) {
    for(size_t n = 0; n < half_length; n+=8){
      //  Grab Twiddle Factors
      __m256d twiddles0 = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
      __m256d twiddles1 = _mm256_load_pd((double*)&local_table[n+2]); // tmp = [r0,i0,r1,i1]
      __m256d twiddles2 = _mm256_load_pd((double*)&local_table[n+4]); // tmp = [r0,i0,r1,i1]
      __m256d twiddles3 = _mm256_load_pd((double*)&local_table[n+6]); // tmp = [r0,i0,r1,i1]
      dual_fft_forward_butterfly(twiddles0, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
      dual_fft_forward_butterfly(twiddles1, (__m256d*)(T+n+2), (__m256d*)(T+n+half_length+2));
      dual_fft_forward_butterfly(twiddles2, (__m256d*)(T+n+4), (__m256d*)(T+n+half_length+4));
      dual_fft_forward_butterfly(twiddles3, (__m256d*)(T+n+6), (__m256d*)(T+n+half_length+6));
    }
  } else {
    // half length may be odd here, so we need to read unaligned single.
    for(size_t n = 0; n < half_length; n++){
      //  Grab Twiddle Factors
      __m128d twiddle = _mm_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
      single_fft_forward_butterfly(twiddle, (__m128d*)(T+n), (__m128d*)(T+n+half_length));
    }
  }

  if(length >= FFT_THRESHOLD_LENGTH) {
    cilk_spawn _fft_forward(T, half_length);
    _fft_forward(T + half_length, half_length);
    cilk_sync;
  } else {
    _fft_forward(T, half_length);
    _fft_forward(T + half_length, half_length);
  }
}

// Runs FFT without cached twiddles
// FIXME: Merge?
void _fft_forward_uncached(complex double *T, size_t length) {
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

  if(table_select(half_length) < twiddle_table_size) {
    if(length >= FFT_THRESHOLD_LENGTH) {
      cilk_spawn _fft_forward(T, half_length);
      _fft_forward(T + half_length, half_length);
      cilk_sync;
    } else {
      _fft_forward(T, half_length);
      _fft_forward(T + half_length, half_length);
    }
  } else {
    if(length >= FFT_THRESHOLD_LENGTH) {
      cilk_spawn _fft_forward_uncached(T, half_length);
      _fft_forward_uncached(T + half_length, half_length);
      cilk_sync;
    } else {
      _fft_forward_uncached(T, half_length);
      _fft_forward_uncached(T + half_length, half_length);
    }
  }
}

void _fft_inverse(complex double *T, size_t length) {
  if (length == 2) {
    dft_2p(T);
    return;
  }

  if(length == 3) {
    dft_3p(T);
    return;
  }

  size_t half_length = length / 2;

  if(length >= FFT_THRESHOLD_LENGTH) {
    cilk_spawn _fft_inverse(T, half_length);
    _fft_inverse(T + half_length, half_length);
    cilk_sync;
  } else {
    _fft_inverse(T, half_length);
    _fft_inverse(T + half_length, half_length);
  }

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[table_select(length)];

  //  Perform FFT reduction into two halves.
  for (size_t n = 0; n < half_length-half_length%2; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddle = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    dual_fft_inverse_butterfly(twiddle, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }
  if(half_length%2 == 1) {
    __m128d twiddle = _mm_load_pd((double*)&local_table[half_length-1]);
    single_fft_inverse_butterfly(twiddle, (__m128d*)(T+half_length-1), (__m128d*)(T+half_length-1+half_length));
  }
}

void _fft_inverse_uncached(complex double *T, size_t length) {
  double omega = 2 * M_PI / length;
  size_t half_length = length / 2;

  if(table_select(half_length) < twiddle_table_size) {
    if(length >= FFT_THRESHOLD_LENGTH) {
      cilk_spawn _fft_inverse(T, half_length);
      _fft_inverse(T + half_length, half_length);
      cilk_sync;
    } else {
      _fft_inverse(T, half_length);
      _fft_inverse(T + half_length, half_length);
    }
  } else {
    if(length >= FFT_THRESHOLD_LENGTH) {
      cilk_spawn _fft_inverse_uncached(T, half_length);
      _fft_inverse_uncached(T + half_length, half_length);
      cilk_sync;
    } else {
      _fft_inverse_uncached(T, half_length);
      _fft_inverse_uncached(T + half_length, half_length);
    }
  }

  //  Perform FFT reduction into two halves.
  for(size_t c = 0; c < half_length; c++) {
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

// External interface
void fft_forward(complex double *T, size_t length) {
  if(table_select(length) < twiddle_table_size) {
    _fft_forward(T, length);
  } else {
    _fft_forward_uncached(T, length);
  }
}

// External interface
void fft_inverse(complex double *T, size_t length) {
  if(table_select(length) < twiddle_table_size) {
    _fft_inverse(T, length);
  } else {
    _fft_inverse_uncached(T, length);
  }
}

