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
#define BAILEYS_THRESHOLD_LENGTH (30) // FIXME: Moar. Just for debugging.

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
    /* FIXME: Implement 5 pt dft */
    // 5*2^k transform
    /*
    if(5*length >= source_length) {
      return 5*length;
    }
    */
    length *= 2;
  }
}

void _fft_forward(complex double *T, size_t length) {
  /*
  if(length == 8) {
    dft_8p(T);
    return;
  }
  */

  /*if(length == 4) {
    dft_4p(T);
    return;
  }*/

  /*
  if(length == 5) {
    dft_5p(T);
    return;
  }
  */

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
  complex double* local_table = twiddle_table[table_select(half_length)];

  for(size_t n = 0; n < half_length-half_length%2; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddles0 = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    dual_fft_forward_butterfly(twiddles0, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }
  if(half_length%2) {
    __m128d twiddle = _mm_load_pd((double*)&local_table[half_length-1]);
    single_fft_forward_butterfly(twiddle, (__m128d*)(T+half_length-1), (__m128d*)(T+half_length-1+half_length));
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

void _strided_fft_forward(complex double *T, size_t length, size_t stride, size_t shift) {
  //if(length == 4) {
  //  dft_4p(T); // FIXME: Stride + Shift
  //  return;
  //}

  size_t content_length = length / stride;

  if(content_length == 2) {
    dft_2p_strided(T, stride, shift);
    return;
  }

  if(content_length == 3) {
    fprintf(stderr, "Strided FFT cannot handle factor 3");
    abort(); // we make sure to not have factor 3 in the strided fft.
    //dft_3p(T); // FIXME: Stride + Shift
    //return;
  }


  size_t half_length = content_length / 2;

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[table_select(half_length)];

  // Dual FFT is tricky, because of the interleaving. Needs another algo (or
  // handle shift and shift+1 at the same time? => Move code to baileys) FIXME
  for(size_t n = 0; n < half_length; n++){
    __m128d twiddle = _mm_load_pd((double*)&local_table[n]);
    single_fft_forward_butterfly(twiddle, (__m128d*)(T+n*stride+shift), (__m128d*)(T+(n+half_length)*stride+shift));
  }

  if(length >= FFT_THRESHOLD_LENGTH) {
    cilk_spawn _strided_fft_forward(T, length/2, stride, shift);
    _strided_fft_forward(T + length/2, length/2, stride, shift);
    cilk_sync;
  } else {
    _strided_fft_forward(T, length/2, stride, shift);
    _strided_fft_forward(T + length/2, length/2, stride, shift);
  }
}

void _baileys_forward(complex double *T, size_t length) {
  size_t a, b;

  b = length;
  a = 1;

  // FIXME: log2, divide by 2, then shift or even faster algo?
  while(a < b) {
    a*=2;
    b/=2;
  }

  for(size_t i=1;i<b;i++) { // FIXME: cilk_for
    _strided_fft_forward(T, length, b, i);
  }


  // Twiddle multiplication
  // FIXME: Could be done in one step together with the strided FFT.
  complex double* local_table = twiddle_table[table_select(length/2)];
  // we can start from 1, because i=0 row has twiddle = 1 (even with bit
  // reverse), so we can skip multiplication
  for(size_t i=0;i<a;i++) {
    size_t istar = bitreverse(i, bitlog2(a));

    // we can start from 1, because j=0 column has twiddle = 1, so we can skip multiplication
    for(size_t j=0;j<b;j++) {
      // Table contains e^(i*2*pi/length*N) for N = [0..length/2]
      // We need value for N = istar*j. We can use %length. If N still > length/2, use complex conj.
      //
      // FIXME: Optimize. Totally shitty sparse access. Either calc from scracth or generate additional tables?
      size_t tblidx = (istar*j)%length;
      int conjug = 0;
      if(tblidx > length/2) {
        tblidx = length - tblidx;
        conjug = 1;
      }
      if(tblidx == length/2) {
        tblidx = 0;
      }
      __m128d twiddle = _mm_load_pd((double*)&local_table[tblidx]);
      if(conjug) {
        twiddle = _mm_xor_pd(twiddle, _mm_set_pd(-0.0, 0.0)); // conjugate.
      }

#ifndef HEAVY_DEBUG
      double real = ((double*)&local_table[tblidx])[0];
      double imag = ((double*)&local_table[tblidx])[1];
      printf("%ld Matrix %ld,%ld (%ld) (which is really %ld,%ld (%ld)), twiddle: %f,%f (CONJ=%d)\n", tblidx, i, j, i*b+j, istar, j, istar*b+j, real, imag, conjug);
      if(tblidx >= length/2) {
        // Sanity check.
        abort();
      }
#endif

      __m128d r0       = _mm_unpacklo_pd(twiddle, twiddle);     //   r = [r0,r0]
      __m128d i0       = _mm_unpackhi_pd(twiddle, twiddle);     //   i = [i0,i0]

      __m128d e0       = _mm_load_pd((double*)&T[i*b+j]);

      // Complex multiplication (result: [er*r0 - ei*i0, er*i0 + ei*r0])
      __m128d a0 = _mm_mul_pd(e0, r0);          //  = [er*r0, ei*r0]
      __m128d b0 = _mm_shuffle_pd(e0, e0, 0x1); //  = [ei, er]
      __m128d c0 = _mm_mul_pd(b0, i0);          //  = [ei*i0, er*i0]
      __m128d d0 = _mm_addsub_pd(a0, c0);       //  = [er*r0 - ei*i0, ei*r0 + er*i0]

      _mm_store_pd((double*)&T[i*b+j], d0);
    }
  }

  for(size_t i=0;i<a;i++) { // FIXME: cilk_for
    _fft_forward(T+i*b, b);
  }

  return;
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
  /*if(length == 4) {
    dft_4p_inv(T);
    return;
  }*/

  if (length == 2) {
    dft_2p(T);
    return;
  }

  if(length == 3) {
    dft_3p_inv(T);
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
  complex double* local_table = twiddle_table[table_select(half_length)];

  //  Perform FFT reduction into two halves.
  for (size_t n = 0; n < half_length-half_length%2; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddle = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    dual_fft_inverse_butterfly(twiddle, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }
  if(half_length%2) {
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
    if(length > BAILEYS_THRESHOLD_LENGTH) {
      _baileys_forward(T, length);
    }
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

