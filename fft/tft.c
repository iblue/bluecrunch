#include <complex.h>
#include <math.h>
#include <x86intrin.h>
#include <stdint.h>
#include <stdio.h> // fprintf
#include <assert.h>
#include "fft.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

// FIXME: Deduplication!!!
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

// FIXME: Deduplication!!!
static inline void dft_2p(complex double* V) {
  __m128d *T = (__m128d*)V;
  __m128d a = T[0];
  __m128d b = T[1];
  T[0] = _mm_add_pd(a,b);
  T[1] = _mm_sub_pd(a,b);
}

static inline int bitlog2(int N) {
  int k = 0;
  while (N >>= 1) {
    k++;
  }
  return k;
}

complex double omega(int i, int N) {
  int k = bitlog2(N);

  complex double* local_table = twiddle_table[k];
  complex double  val         = local_table[i];

  return val;
}

void tft_inverse1(complex double *T, size_t head, size_t tail, size_t last, size_t s) {
  size_t left_middle  = (last - head)/2 + head;
  size_t right_middle = left_middle + 1;

  if(head > tail) {
    return;
  } else if(tail >= left_middle) {
    // Recursion end:
    if(last-head==1) {
      if(tail-head==0) {
        T[head] = 2*T[head] - T[last];
        return;
      } else {
        fprintf(stderr, "Not implemented\n");
        return;
      }
    }

    // Push up the self-contained region T[head] to T[left_middle]
    fft_inverse(T+head, bitlog2(left_middle - head + 1), 1);

    // Push down T[tail+1] .. T[last] from T[tail+1-left_middle] .. T[left_middle]
    for(size_t i=tail+1;i<=last;i++) {
      complex double c = T[head+i-left_middle-1];
      complex double b = T[i];
      T[i]  = c - b;
      T[i] *= omega(i-left_middle-1, last-head+1);
    }

    tft_inverse1(T, right_middle, tail, last, s+1);

    size_t m_s = bitlog2(last-head+1);
    size_t half_length = 1 << (m_s - 1);
    complex double* local_table = twiddle_table[m_s];

    // Push up (in pairs) (T[head], T[head+m_s]) ... (T[left_middle], T[left_middle+m_s])
    if(half_length == 1) {
      dft_2p(T+head);
    } else for(size_t n=0; n < half_length; n+=2) {
      __m256d twiddle = _mm256_load_pd((double*)&local_table[n]);
      fft_inverse_butterfly(twiddle, (__m256d*)(T+head+n), (__m256d*)(T+head+n+half_length));
    }
  } else if(tail < left_middle) {
    fprintf(stderr, "Not implemented\n");
    abort();
  }
}

void tft_inverse(complex double *T, size_t len) {
  size_t n = 2;
  size_t l = len;

  assert(len >= 2);

  while(len /= 2) {
    n *= 2;
  }

  tft_inverse1(T, 0, l-1, n-1, 1);
}
