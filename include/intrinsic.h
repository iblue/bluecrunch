#include <x86intrin.h>

// Returns ceil(log2(N))
static inline __attribute__((always_inline)) int bitlog2(int N) {
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

static inline __attribute__((always_inline)) size_t bitreverse(size_t index, size_t bits) {
  static const unsigned char tbl[] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
    0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
    0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
    0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
    0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
    0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
    0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
    0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
    0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
    0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
    0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
    0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
    0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
    0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
    0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
    0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
  };

  size_t returns=0;

  unsigned char *in  = (unsigned char *)&(index);
  unsigned char *out = (unsigned char *)&(returns);

  // TODO: Performance? Stop after n bits, because we can be sure that index <
  // 2^bits, so the rest must be 0.
  for(int i=0;i<8;i++) {
    out[7-i] = tbl[in[i]];
  }

  returns >>= 64-bits;

  return returns;
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
static inline __attribute__((always_inline)) void dual_fft_forward_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
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

static inline __attribute__((always_inline)) void single_fft_forward_butterfly(__m128d twiddle, __m128d* T0, __m128d* T1) {
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
static inline __attribute__((always_inline)) void dual_fft_inverse_butterfly(__m256d twiddles, __m256d* T0, __m256d* T1) {
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

static inline __attribute__((always_inline)) void single_fft_inverse_butterfly(__m128d twiddle, __m128d* T0, __m128d* T1) {
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

static inline __attribute__((always_inline)) void dft_2p(complex double* V) {
  __m128d *T = (__m128d*)V;
  __m128d a = T[0];
  __m128d b = T[1];
  T[0] = _mm_add_pd(a,b);
  T[1] = _mm_sub_pd(a,b);
}

static inline __attribute__((always_inline)) void dft_2p_strided(complex double* V, size_t stride, size_t shift) {
  __m128d *T = (__m128d*)V;
  __m128d a = T[shift];
  __m128d b = T[stride+shift];
  T[shift] = _mm_add_pd(a,b);
  T[stride+shift] = _mm_sub_pd(a,b);
}

#define OMEGA_3_R (0.5)
#define OMEGA_3_I (0.8660254037844386467637231707529361834714026269051903140)
/*
 * So we want to calculate:
 *
 * complex double omega_1_3 = OMEGA_3_R + I*OMEGA_3_I;
 * complex double omega_2_3 = OMEGA_3_R - I*OMEGA_3_I;
 *
 * complex double a = V[0];
 * complex double b = V[1];
 * complex double c = V[2];
 *
 * V[0] = a + b + c;
 * V[1] = a + omega_1_3 * b + omega_2_3 * c;
 * V[2] = a + omega_2_3 * b + omega_1_3 * c;
 *
 * Complex multiplication:
 * let x = a + bi = [a, b]
 * let y = c + di = [c, d]
 *
 * x*y   = [ac-bd, (a+b)(c+d) - ac - bd] OR
 *         [ac-bd, ad+bc]
 *
 * With:
 * v[0] = [ar, ai] +            [br, bi] +            [cr, ci]
 * v[1] = [ar, ai] + [-or,  oi]*[br, bi] + [-or, -oi]*[cr, ci]
 * v[2] = [ar, ai] + [-or, -oi]*[br, bi] + [-or,  oi]*[cr, ci]
 *
 * We get:
 * v[0] = [ar, ai] +            [br, bi] +            [cr, ci]
 * v[1] = [ar, ai] + [-or*br-oi*bi, (br+bi)(-or+oi)+or*br-oi*bi]
 *                 + [-or*cr+oi*ci, (cr+ci)(or+oi) +or*cr+oi*ci]
 * v[2] = [ar, ai] + [-or*br+oi*bi, (br+bi)(-or-oi)+or*br+oi*bi]
 *                 + [-or*cr-oi*ci, (cr+ci)(or-oi) +or*cr-oi*ci]
 *
 * OR
 *
 * v[0] = [ar, ai] +                     [br, bi] +            [cr, ci]
 * v[1] = [ar, ai] + [-or*br-oi*bi, -or*bi+oi*br] + [-or*cr+oi*ci, -or*ci-oi*cr]
 * v[2] = [ar, ai] + [-or*br+oi*bi, -or*bi-oi*br] + [-or*cr-oi*ci, -or*ci+oi*cr]
 *
 * Operations:
 * - addsub
 * - shuffle
 * - mul
 *
 * or  = [-or, -or]
 * oi  = [oi, oi]
 * b   = [br, bi]
 * c   = [cr, ci]
 *
 * // calculate B part
 * n0 = or*b     // mulpd:   n0 = [-or*br, -or*bi]
 * n1 = shuff b  // shuffle  n1 = [bi, br]
 * n2 = oi*n1    // mulpd    n2 = [oi*bi, oi*br]
 * v2 = n0+-n2   // addsub   v2 = [-or*br+oi*bi, -or*bi-oi*br]
 * n3 = -n2      // xor      n3 = [-oi*bi, -oi*br]
 * v1 = n0+-n3   // addsub   v1 = [-or*br-oi*bi, -or*bi+oi*br]
 *
 * // calculate C part
 * m0 = or*c     // mulpd:   m0 = [-or*cr, -or*ci]
 * m1 = shuff c  // shuffle  m1 = [ci, cr]
 * m2 = oi*m1    // mulpd    m2 = [oi*ci, oi*cr]
 * w1 = m0+-m2   // addsub   w2 = [-or*cr+oi*ci, -or*ci-oi*cr]
 * m3 = -m2      // xor      m3 = [-oi*ci, -oi*cr]
 * w2 = m0+-m3   // addsub   w1 = [-or*cr-oi*ci, -or*ci+oi*cr]
 *
 * // add parts
 * V0 = a + b + c
 * V1 = a + v1 + w1
 * V2 = a + v2 + w2
 */
static inline __attribute__((always_inline)) void dft_3p(complex double* V) {
#ifdef HEAVY_DEBUG
  complex double omega_1_3 = -OMEGA_3_R + I*OMEGA_3_I;
  complex double omega_2_3 = -OMEGA_3_R - I*OMEGA_3_I;

  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];

  complex double x = a + b + c;
  printf("expected: V[0] = [%f, %f]\n", creal(x), cimag(x));
  complex double y = a + omega_1_3 * b + omega_2_3 * c;
  printf("expected: V[1] = [%f, %f]\n", creal(y), cimag(y));
  complex double z = a + omega_2_3 * b + omega_1_3 * c;
  printf("expected: V[2] = [%f, %f]\n", creal(z), cimag(z));

  complex double M11 = omega_1_3 * b;
  complex double M12 = omega_2_3 * c;
  complex double M21 = omega_2_3 * b;
  complex double M22 = omega_1_3 * c;
  printf("expected: M11 = [%f, %f]\n", creal(M11), cimag(M11));
  printf("expected: M12 = [%f, %f]\n", creal(M12), cimag(M12));
  printf("expected: M21 = [%f, %f]\n", creal(M21), cimag(M21));
  printf("expected: M22 = [%f, %f]\n", creal(M22), cimag(M22));
#endif

  __m128d orm = _mm_set1_pd(-OMEGA_3_R);
  __m128d oim = _mm_set1_pd(OMEGA_3_I);

  __m128d am = _mm_load_pd((double*)V);
  __m128d bm = _mm_load_pd((double*)(V+1));
  __m128d cm = _mm_load_pd((double*)(V+2));

  __m128d n0 = _mm_mul_pd(orm, bm);
  __m128d n1 = _mm_shuffle_pd(bm, bm, 0x1);
  __m128d n2 = _mm_mul_pd(oim, n1);
  __m128d m11 = _mm_addsub_pd(n0, n2);
  __m128d n3 = _mm_xor_pd(n2, _mm_set1_pd(-0.0));
  __m128d m21 = _mm_addsub_pd(n0, n3);

  __m128d m0 = _mm_mul_pd(orm, cm);
  __m128d m1 = _mm_shuffle_pd(cm, cm, 0x1);
  __m128d m2 = _mm_mul_pd(oim, m1);
  __m128d m22 = _mm_addsub_pd(m0, m2);
  __m128d m3 = _mm_xor_pd(m2, _mm_set1_pd(-0.0));
  __m128d m12 = _mm_addsub_pd(m0, m3);

  __m128d V0 = _mm_add_pd(am, bm);
          V0 = _mm_add_pd(V0, cm);
  __m128d V1 = _mm_add_pd(am, m11);
          V1 = _mm_add_pd(V1, m12);
  __m128d V2 = _mm_add_pd(am, m21);
          V2 = _mm_add_pd(V2, m22);

  // FIXME: only real inputs? then V2 is complex conjugate of V1 (first round)

  _mm_store_pd((double*)V,     V0);
  _mm_store_pd((double*)(V+1), V1);
  _mm_store_pd((double*)(V+2), V2);
}

static inline __attribute__((always_inline)) void dft_3p_inv(complex double* V) {
  __m128d orm = _mm_set1_pd(-OMEGA_3_R);
  __m128d oim = _mm_set1_pd(-OMEGA_3_I);

  __m128d am = _mm_load_pd((double*)V);
  __m128d bm = _mm_load_pd((double*)(V+1));
  __m128d cm = _mm_load_pd((double*)(V+2));

  __m128d n0 = _mm_mul_pd(orm, bm);
  __m128d n1 = _mm_shuffle_pd(bm, bm, 0x1);
  __m128d n2 = _mm_mul_pd(oim, n1);
  __m128d m11 = _mm_addsub_pd(n0, n2);
  __m128d n3 = _mm_xor_pd(n2, _mm_set1_pd(-0.0));
  __m128d m21 = _mm_addsub_pd(n0, n3);

  __m128d m0 = _mm_mul_pd(orm, cm);
  __m128d m1 = _mm_shuffle_pd(cm, cm, 0x1);
  __m128d m2 = _mm_mul_pd(oim, m1);
  __m128d m22 = _mm_addsub_pd(m0, m2);
  __m128d m3 = _mm_xor_pd(m2, _mm_set1_pd(-0.0));
  __m128d m12 = _mm_addsub_pd(m0, m3);

  __m128d V0 = _mm_add_pd(am, bm);
          V0 = _mm_add_pd(V0, cm);
  __m128d V1 = _mm_add_pd(am, m11);
          V1 = _mm_add_pd(V1, m12);
  __m128d V2 = _mm_add_pd(am, m21);
          V2 = _mm_add_pd(V2, m22);

  _mm_store_pd((double*)V,     V0);
  _mm_store_pd((double*)(V+1), V1);
  _mm_store_pd((double*)(V+2), V2);
}

/*
  * V[0] = a +   b + c +   d; // [ar+br+cr+dr, ai+bi+ci+di]
  * V[1] = a -   b + c -   d; // [ar-br+cr-dr, ai-bi+ci-di]
  * V[2] = a + I*b - c - I*d; // [ar-bi-cr+di, ai+br-ci-dr]
  * V[3] = a - I*b - c + I*d; // [ar+bi-cr-di, ai-br-ci+dr]
  */
static inline __attribute__((always_inline)) void dft_4p(complex double* V) {
#ifdef HEAVY_DEBUG
  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];
  complex double d = V[3];

  complex double x = a + b + c + d;
  printf("expected: V[0] = [%f, %f]\n", creal(x), cimag(x));
  complex double y = a - b + c - d;
  printf("expected: V[1] = [%f, %f]\n", creal(y), cimag(y));
  complex double z = a + I*b - c - I*d;
  printf("expected: V[2] = [%f, %f]\n", creal(z), cimag(z));
  complex double q = a - I*b - c + I*d;
  printf("expected: V[3] = [%f, %f]\n", creal(q), cimag(q));
#endif

  // load
  __m128d am = _mm_load_pd((double*)V);
  __m128d bm = _mm_load_pd((double*)(V+1));
  __m128d cm = _mm_load_pd((double*)(V+2));
  __m128d dm = _mm_load_pd((double*)(V+3));

  // for negating
  __m128d gg = _mm_set1_pd(-0.0);

  // v1, v2
  __m128d n0 = _mm_add_pd(am, cm); // n0 = [ar+cr, ai+ci]
  __m128d n1 = _mm_add_pd(bm, dm); // n1 = [br+dr, bi+di]
  __m128d v0 = _mm_add_pd(n0, n1); // v0 = [ar+br+cr+dr, ai+bi+ci+di]
  __m128d v1 = _mm_sub_pd(n0, n1); // v1 = [ar-br+cr-dr, ai-bi+ci-di]

  // v3, v4
  __m128d n2 = _mm_sub_pd(am, cm); // n2 = [ar-cr, ai-ci];
  __m128d n3 = _mm_shuffle_pd(bm, bm, 0x1); // n3 = [bi, br]
  __m128d z0 = _mm_addsub_pd(n2, n3); // z0 = [ar+bi-cr, ai-br-ci] (-> V3)
  __m128d n4 = _mm_shuffle_pd(dm, dm, 0x1); // n4 = [di, dr]
  __m128d z1 = _mm_addsub_pd(n2, n4); // z1 = [ar-cr+di, ai-ci-dr] (-> V2)
  __m128d q0 = _mm_xor_pd(n3, gg);    // q0 = [-bi, -br]
  __m128d v3 = _mm_addsub_pd(z1, q0); // v3 = [ar-bi-cr+di, ai+br-ci-dr]
  __m128d q1 = _mm_xor_pd(n4, gg);    // q1 = [-di, -dr]
  __m128d v2 = _mm_addsub_pd(z0, q1); // v2 = [ar+bi-cr-di, ai-br-ci+dr]

  // store
  _mm_store_pd((double*)V,     v0);
  _mm_store_pd((double*)(V+1), v1);
  _mm_store_pd((double*)(V+2), v2);
  _mm_store_pd((double*)(V+3), v3);
}

/*
 * V[0] = a + b +   c +   d; // [ar+br+cr+dr, ai+bi+ci+di]
 * V[1] = a - b - I*c + I*d; // [ar-br+ci-di, ai-bi-cr+dr]
 * V[2] = a + b -   c -   d; // [ar+br-cr-dr, ai+bi-ci-di]
 * V[3] = a - b + I*c - I*d; // [ar-br-ci+di, ai-bi+cr-dr]
 *
 */
static inline __attribute__((always_inline)) void dft_4p_inv(complex double* V) {
#ifdef HEAVY_DEBUG
  complex double a = V[0];
  complex double b = V[1];
  complex double c = V[2];
  complex double d = V[3];

  complex double x = a + b +   c +   d; // [ar+br+cr+dr, ai+bi+ci+di]
  complex double y = a - b - I*c + I*d; // [ar-br+ci-di, ai-bi-cr+dr]
  complex double z = a + b -   c -   d; // [ar+br-cr-dr, ai+bi-ci-di]
  complex double q = a - b + I*c - I*d; // [ar-br-ci+di, ai-bi+cr-dr]
  printf("expected: V[0] = [%f, %f]\n", creal(x), cimag(x));
  printf("expected: V[1] = [%f, %f]\n", creal(y), cimag(y));
  printf("expected: V[2] = [%f, %f]\n", creal(z), cimag(z));
  printf("expected: V[3] = [%f, %f]\n", creal(q), cimag(q));
#endif

  // load
  __m128d am = _mm_load_pd((double*)V);
  __m128d bm = _mm_load_pd((double*)(V+1));
  __m128d cm = _mm_load_pd((double*)(V+2));
  __m128d dm = _mm_load_pd((double*)(V+3));

  // for negating
  __m128d gg = _mm_set1_pd(-0.0);

  // v0, v2
  __m128d n0 = _mm_add_pd(am, bm);  // n0 = [ar+br, ai+bi]
  __m128d n1 = _mm_add_pd(cm, dm);  // n1 = [cr+dr, ci+di]
  __m128d v0 = _mm_add_pd(n0, n1);  // v0 = [ar+br+cr+dr, ai+bi+ci+di]
  __m128d v2 = _mm_sub_pd(n0, n1);  // v1 = [ar+br-cr-dr, ai+bi-cr-dr]

  // v1, v3
  __m128d ff = _mm_sub_pd(am, bm);          // ff = [ar-br, ai-bi]
  __m128d z0 = _mm_sub_pd(cm, dm);          // z0 = [cr-dr, ci-di]*/
  __m128d q0 = _mm_shuffle_pd(z0, z0, 0x1); // q0 = [ci-di, cr-dr]
  __m128d v3 = _mm_addsub_pd(ff, q0);       // v3 = [ar-br-ci+di, ai-bi+cr-dr]
  __m128d q1 = _mm_xor_pd(q0, gg);          // q1 = [-ci+di, -cr+dr]
  __m128d v1 = _mm_addsub_pd(ff, q1);       // v1 = [ar-br+ci-di, ai-bi-cr+dr]

  // sequential access
  _mm_store_pd((double*)V, v0);
  _mm_store_pd((double*)(V+1), v1);
  _mm_store_pd((double*)(V+2), v2);
  _mm_store_pd((double*)(V+3), v3);

#ifdef HEAVY_DEBUG
  if(V[0] != x || V[1] != y || V[2] != z || V[3] != q) {
    printf("DFT4 inv error\n");
    abort();
  }
#endif
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

