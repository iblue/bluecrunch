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

