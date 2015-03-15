#include <x86intrin.h>
#include <stdint.h>
#include "fft.h"

// Inline pointwise multiply of complex numbers
// Input: 2 vectors of complex numbers
// Output: Sets first vector to product
void fft_pointwise(__m128d *a, __m128d *b, int k) {

  // FFT length
  size_t len = 1 << k;

  for (size_t i = 0; i < len; i++) {
    // TODO: AVX
    // Complex multiplication in SSE3
    __m128d a0 = a[i];
    __m128d b0 = b[i];
    __m128d c0 = _mm_mul_pd(a0, _mm_unpacklo_pd(b0,b0));
    __m128d d0 = _mm_mul_pd(_mm_shuffle_pd(a0, a0, 1), _mm_unpackhi_pd(b0, b0));

    a[i] = _mm_addsub_pd(c0, d0);
  }
}

// Inline pointwise multiply of complex numbers
// Input: 2 vectors of complex numbers
// Output: Sets first vector to product
void tft_pointwise(__m128d *a, __m128d *b, size_t len) {
  for (size_t i = 0; i < len; i++) {
    // TODO: AVX
    // Complex multiplication in SSE3
    __m128d a0 = a[i];
    __m128d b0 = b[i];
    __m128d c0 = _mm_mul_pd(a0, _mm_unpacklo_pd(b0,b0));
    __m128d d0 = _mm_mul_pd(_mm_shuffle_pd(a0, a0, 1), _mm_unpackhi_pd(b0, b0));

    a[i] = _mm_addsub_pd(c0, d0);
  }
}
