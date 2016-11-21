#include <x86intrin.h>

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

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

// Runs a forward butterfly on 2 elements simultaneously
// Arguments:
// - twiddles: Vector of twiddles [r0, i0, r1, i1]
// - T0:       Pointer to 2 complex doubles (first half)
// - T1:       Pointer to 2 complex doubles (second half)
//
//  Perform FFT reduction into two halves.
// Input: a[0], a[1], b[0], b[1], W[0], W[1]
//        ----------  ----------  ----------
//            |           |           |
//      real & imag T[c]  |    real & imag root
//               real & imag T[c+N/2]
//
// Output:
//  a[0] <- a[0] + b[0]
//  a[1] <- a[1] + b[1]
//  b[0] <- W[0]*(a[0] - b[0]) - W[1]*(a[1] - b[1])
//  b[1] <- W[0]*(a[1] - b[1]) - W[1]*(a[0] - b[0])
static inline void dual_fft_forward_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
  __m256d r   = _mm256_unpacklo_pd(twiddles, twiddles);  //   r = [r0,r0,r1,r1]
  __m256d i   = _mm256_unpackhi_pd(twiddles, twiddles);  //   i = [i0,i0,i1,i1]

  //  Grab elements
  __m256d a = _mm256_load_pd((double*)T0);  // a = [a0r, a0i, a1r, a1i]
  __m256d b = _mm256_loadu_pd((double*)T1);  // b = [b0r, b0i, b1r, b1i] (may be misaligned)

  //  Perform butterfly
  __m256d c = _mm256_add_pd(a,b); // c = a + b
  __m256d d = _mm256_sub_pd(a,b); // d = a - b

  _mm256_store_pd((double*)T0, c); // T[c] <- c

  //  Multiply by twiddle factor.
  //  FIXME: Documentation. Tomorrow will have forgotten how this works!
  c = _mm256_mul_pd(d,r);                             // c = (a-b)*omega_r
  d = _mm256_mul_pd(_mm256_shuffle_pd(d,d,0x5), i);   // shuffle -> switch real and imag
  c = _mm256_addsub_pd(c,d);

  // (may be misaligned)
  _mm256_storeu_pd((double*)T1, c); // T[c+N/2] <- c
}

static inline void single_fft_forward_butterfly(__m128d twiddle, __m128d* T0, __m128d* T1) {
  __m128d r   = _mm_unpacklo_pd(twiddle, twiddle);  //   r = [r0,r0]
  __m128d i   = _mm_unpackhi_pd(twiddle, twiddle);  //   i = [i0,i0]

  //  Grab elements
  __m128d a = _mm_load_pd((double*)T0);  // a = [a0r, a0i]
  __m128d b = _mm_load_pd((double*)T1);  // b = [b0r, b0i]

  //  Perform butterfly
  __m128d c = _mm_add_pd(a,b); // c = a + b
  __m128d d = _mm_sub_pd(a,b); // d = a - b

  _mm_store_pd((double*)T0, c); // T[c] <- c

  //  Multiply by twiddle factor.
  //  FIXME: Documentation. Tomorrow will have forgotten how this works!
  c = _mm_mul_pd(d,r);                             // c = (a-b)*omega_r
  d = _mm_mul_pd(_mm_shuffle_pd(d,d,0x1), i);   // shuffle -> switch real and imag
  c = _mm_addsub_pd(c,d);

  _mm_store_pd((double*)T1, c); // T[c+N/2] <- c
}

// Runs an inverse butterfly on 2 elements simultaneously
// Arguments:
// - twiddles: Vector of twiddles [r0, i0, r1, i1]
// - T0:       Pointer to 2 complex doubles (first half)
// - T1:       Pointer to 2 complex doubles (second half)
static inline void dual_fft_inverse_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
  __m256d r   = _mm256_unpacklo_pd(twiddles, twiddles);     //   r = [r0,r0,r1,r1]
  __m256d i   = _mm256_unpackhi_pd(twiddles, twiddles);     //   i = [i0,i0,i1,i1]
  i = _mm256_xor_pd(i,_mm256_set1_pd(-0.0));              //   i = -i

  //  Grab elements
  __m256d a = _mm256_load_pd((double*)T0);  // a = [a0r, a0i, a1r, a1i]
  __m256d b = _mm256_loadu_pd((double*)T1);  // b = [b0r, b0i, b1r, b1i] (maybe misaligned)

  //  Perform butterfly
  __m256d c, d;

  //  Multiply by twiddle factor.
  c = _mm256_mul_pd(b,r);
  d = _mm256_mul_pd(_mm256_shuffle_pd(b,b,0x5),i);
  c = _mm256_addsub_pd(c,d);

  b = _mm256_add_pd(a,c);
  d = _mm256_sub_pd(a,c);

  _mm256_store_pd((double*)T0, b);
  _mm256_storeu_pd((double*)T1, d); // (maybe misaligned)
}

static inline void single_fft_inverse_butterfly(__m128d twiddle, __m128d* T0, __m128d* T1) {
  __m128d r   = _mm_unpacklo_pd(twiddle, twiddle);     //   r = [r0,r0]
  __m128d i   = _mm_unpackhi_pd(twiddle, twiddle);     //   i = [i0,i0]
  i = _mm_xor_pd(i,_mm_set1_pd(-0.0));              //   i = -i

  //  Grab elements
  __m128d a = _mm_load_pd((double*)T0);  // a = [a0r, a0i]
  __m128d b = _mm_load_pd((double*)T1);  // b = [b0r, b0i]

  //  Perform butterfly
  __m128d c, d;

  //  Multiply by twiddle factor.
  c = _mm_mul_pd(b,r);
  d = _mm_mul_pd(_mm_shuffle_pd(b,b,0x1),i);
  c = _mm_addsub_pd(c,d);

  b = _mm_add_pd(a,c);
  d = _mm_sub_pd(a,c);

  _mm_store_pd((double*)T0, b);
  _mm_store_pd((double*)T1, d);
}

static inline void dft_2p(complex double* V) {
  __m128d *T = (__m128d*)V;
  __m128d a = T[0];
  __m128d b = T[1];
  T[0] = _mm_add_pd(a,b);
  T[1] = _mm_sub_pd(a,b);
}

// FIXME: We surely can optimize this shit with some magic.
static inline void dft_3p(complex double* V) {
  complex double omega_1_3 = -0.5 + 0.8660254037844386467637231707529361834714026269051903140*I;
  complex double omega_2_3 = -0.5 - 0.8660254037844386467637231707529361834714026269051903140*I;

  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];

  V[0] = a + b + c;
  V[1] = a + omega_1_3 * b + omega_2_3 * c;
  V[2] = a + omega_2_3 * b + omega_1_3 * c;
}

// FIXME: We surely can optimize this shit with some magic.
static inline void dft_3p_inv(complex double* V) {
  complex double omega_1_3 = -0.5 + 0.8660254037844386467637231707529361834714026269051903140*I;
  complex double omega_2_3 = -0.5 - 0.8660254037844386467637231707529361834714026269051903140*I;

  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];

  V[0] = a + b + c;
  V[1] = a + omega_2_3 * b + omega_1_3 * c;
  V[2] = a + omega_1_3 * b + omega_2_3 * c;
}

static inline void dft_4p(complex double* V) {
  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];
  complex double d = V[3];

  V[0] = a +   b + c +   d;
  V[1] = a -   b + c -   d;
  V[2] = a + I*b - c - I*d;
  V[3] = a - I*b - c + I*d;
}

// F_8 = B_8 (I_2 (x) F_4)
/*
static inline void dft_8p(complex double* V) {
  double invsqrt_2 = 0.707106781186547524400844362104849039284835937688474036588;

  complex double B[] = {
    1.0,
    invsqrt_2 + invsqrt_2*I,
    I,
    -invsqrt_2 + invsqrt_2*I
  };

  dft_4p(V);
  dft_4p(V+4);

  complex double a = V[0] + B[0]*V[4];
  complex double b = V[1] + B[1]*V[5];
  complex double c = V[2] + B[2]*V[6];
  complex double d = V[3] + B[3]*V[7];
  complex double e = V[0] - B[0]*V[4];
  complex double f = V[1] - B[1]*V[5];
  complex double g = V[2] - B[2]*V[6];
  complex double h = V[3] - B[3]*V[7];

  V[0] = a;
  V[4] = b;
  V[2] = c;
  V[6] = d;
  V[1] = e;
  V[5] = f;
  V[3] = g;
  V[7] = h;
}
*/

// 
// [ 1   1   ] [ 1  1     ] [1      ] [V0]
// [   1   i ] [ 1 -1     ] [    1  ] [V1]
// [ 1   -1  ] [      1  1] [  1    ] [V2]
// [   1   -i] [      1 -1] [      1] [V3]
//
// x' = P_4^T (I_1 (x) B_4) (I_2 (x) B_2) P_4^T x
/*
static inline void dft_4p(complex double* V) {
  // u = PI^T_4 v
  complex double u0 = V[0];
  complex double u1 = V[2];
  complex double u2 = V[1];
  complex double u3 = V[3];

  // w = (I_2 (x) B_2) u
  complex double w0 = u0+u1;
  complex double w1 = u0-u1;
  complex double w2 = u2+u3;
  complex double w3 = u2-u3;

  // v' = PI^T_4 (I_1 (x) B_4) u
  V[0] = w0 + w2;     // = a + b + c + d
  V[2] = w1 + I*w3;   // = u0-u1 + I*(u2-u3) = a+I*b-c-I*d
  V[1] = w0 - w2;     // = a-b+d-c
  V[3] = w1 - I*w3;   // = a-I*b-c+I*d
}
*/

static inline void dft_4p_inv(complex double* V) {
  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];
  complex double d = V[3];

  V[0] = a + b +   c +   d;
  V[1] = a - b - I*c + I*d;
  V[2] = a + b -   c -   d;
  V[3] = a - b + I*c - I*d;
}
