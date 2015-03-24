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
static inline void fft_forward_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
  __m256d r   = _mm256_unpacklo_pd(twiddles, twiddles);  //   r = [r0,r0,r1,r1]
  __m256d i   = _mm256_unpackhi_pd(twiddles, twiddles);  //   i = [i0,i0,i1,i1]

  //  Grab elements
  __m256d a = _mm256_load_pd((double*)T0);  // a = [a0r, a0i, a1r, a1i]
  __m256d b = _mm256_load_pd((double*)T1);  // b = [b0r, b0i, b1r, b1i]

  //  Perform butterfly
  __m256d c = _mm256_add_pd(a,b); // c = a + b
  __m256d d = _mm256_sub_pd(a,b); // d = a - b

  _mm256_store_pd((double*)T0, c); // T[c] <- c

  //  Multiply by twiddle factor.
  //  FIXME: Documentation. Tomorrow will have forgotten how this works!
  c = _mm256_mul_pd(d,r);                             // c = (a-b)*omega_r
  d = _mm256_mul_pd(_mm256_shuffle_pd(d,d,0x5), i);   // shuffle -> switch real and imag
  c = _mm256_addsub_pd(c,d);

  _mm256_store_pd((double*)T1, c); // T[c+N/2] <- c
}

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

// Returns ceil(log2(N))
static inline int bitlog2(int N) {
  // FIXME: Branchless code?

  // Looks like there is no intrinsic for the x86 bsr instruction, but this
  // will roughly generate the following code with optimization turned on:
  // popcnt %esi, %eax
  // bsr    %esi, %ecx
  // cmp    $0x1, %eax
  // ...
  // je     label
  // add    $0x1, %ecx
  // label:
  //
  if(__builtin_popcount(N) == 1) {
    return __builtin_clz(N) ^ 0x1f;
  } else {
    return (__builtin_clz(N) ^ 0x1f) + 1;
  }
}

complex double omega(int i, int N) {
  int k = bitlog2(N);

  complex double* local_table = twiddle_table[k];
  complex double  val         = local_table[i];

  return val;
}

void tft_forward(complex double *T, size_t length) {
  size_t k = bitlog2(length);
  size_t full_length = 1 << k;

  //size_t left_length = full_length / 2;
  //size_t right_length = full_length - left_length;

  // (Bit reversed) 2-point DFT
  if(k==1) {
    dft_2p(T);
    return;
  }

  size_t half_length = full_length / 2;

  //  Get local twiddle table.
  complex double* local_table = twiddle_table[k];

  for (size_t n = 0; n < half_length; n+=2){
    //  Grab Twiddle Factors
    __m256d twiddles = _mm256_load_pd((double*)&local_table[n]); // tmp = [r0,i0,r1,i1]
    fft_forward_butterfly(twiddles, (__m256d*)(T+n), (__m256d*)(T+n+half_length));
  }

  //  No more threads.
  fft_forward(T, k - 1, 1);
  fft_forward(T + half_length, k - 1, 1);
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
    // Push down T[tail+1] to T[left_middle]
    assert(tail+1 <= left_middle); // FIXME?
    for(size_t i=tail+1;i>=left_middle;i--) {
      complex double a = T[i];
      complex double b = T[i+last-left_middle];
      T[i]  = (a+b)/2;
    }

    tft_inverse1(T, head, tail, left_middle, s+1);

    // Push up T[head] to T[left_middle]
    for(size_t i=head;i<=left_middle;i++) {
      complex double b = T[i+last-left_middle];
      T[i] *= 2;
      T[i] -= b;
    }
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
