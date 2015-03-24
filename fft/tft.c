#include <complex.h>
#include <math.h>
#include <x86intrin.h>
#include <stdint.h>
#include <stdio.h> // fprintf
#include <assert.h>
#include "fft.h"
#include "intrinsic.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

complex double omega(int i, int N) {
  int k = bitlog2(N);

  complex double* local_table = twiddle_table[k];
  complex double  val         = local_table[i];

  return val;
}

void tft_forward(complex double *T, size_t length, int k) {
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
      complex double b = 0;

      if(s > 1) {
        b = T[i];
      }

      complex double c = T[head+i-left_middle-1];
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

void tft_inverse(complex double *T, size_t len, int k) {
  size_t n = 1 << k;
  size_t l = len;

  assert(len >= 2);

  // If 2^k normal sized FFT, use FFT algoritmn, because it's faster.
  if(n == l) {
    fft_inverse(T, k, 1);
    return;
  }

  tft_inverse1(T, 0, l-1, n-1, 1);
}
