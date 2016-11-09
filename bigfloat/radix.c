#include <malloc.h>
#include <string.h> // memcpy
#include <stdlib.h> // abort
#include "math.h"
#include "bigfloat.h"

#define OLDBASE (4294967296)
#define NEWBASE (100000000)
#define LOG (log(OLDBASE)/log(NEWBASE))
#define KT (3500)

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

// Calculates approximate number of digits in base 100000000 of a number.
size_t bigfloat_radix_decimals(bigfloat_t target) {
  double approx = 0.0;

  if(target->len >= 3) {
    approx += (double)target->coef[target->len-3];
  }

  if(target->len >= 2) {
    approx += (double)target->coef[target->len-2] * (double)OLDBASE;
  }

  approx += (double)target->coef[target->len-1] * (double)OLDBASE * (double)OLDBASE;

  double digits = log(approx)/log(NEWBASE);
  double additional = ((double)target->len-3.0)*log(OLDBASE)/log(NEWBASE);

  double epsilon = 1e-10;

  return ceil(digits + additional + epsilon);
}

void convert_trunc(bigfloat_t s, const bigfloat_t y0, size_t k, size_t n) {
  double alpha = 1./LOG;

  bigfloat_t t;
  bigfloat_new(t);
  bigfloat_copy(t, y0);
  bigfloat_realloc(t, t->len+1); // To make multiplication succeed

  // FIXME: Performance: realloc is not needed
  // Same probably goes for copy. Move into one new function

  // Save original t ptr for deallocation;
  uint32_t *coef_ptr = t->coef;

  for(size_t i=1;i<=k;i++) {
    //printf("i = %ld\n", i);

    // Multiplication stage: t[i] = b*y[i-1]
    bigfloat_mului(t, NEWBASE);
    //bigfloat_print("t", t);

    // n_{i-1} = n - floor(i*a)
    size_t ni_1 = n - floor(alpha*(double)(i-1));
    size_t ni   = n - floor(alpha*(double)i);
    //printf("ni_1, ni = %ld, %ld\n", ni_1, ni);

    // This is the digit: s[k-i] = t[i] / (2^32)^n_{i-1};
    s->coef[k-i] = t->coef[ni_1];
    //printf("-> %d\n", s->coef[k-i]);

    // Truncate stage: z[i] = t[i] mod 2^{n_{i-1}}
    /*for(size_t j=ni+1;j<t->len;j++) {
      t->coef[j] = 0;
    }*/
    t->len = ni_1;

    // Truncate other side: y[i] = zi bdiv 2^(n_(i-1) - n_i)
    t->coef += (ni_1 - ni);
    t->len  -= (ni_1 - ni);

    // Add a zero to fix multiplication in next iter
    t->len  += 1;
    t->coef[t->len-1] = 0;
  }

  s->len = k+1;

  // Restore ptr and deallocate
  t->coef = coef_ptr;
  bigfloat_free(t);
}

void convert_rec(bigfloat_t s, size_t k, const bigfloat_t y, size_t n, size_t g) {
  //if(k <= KT) {
    convert_trunc(s, y, k, n);
  //} else {
    // TODO
  //}
}

// Converts from OLDBASE to NEWBASE
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  size_t n, k, g;
  if(-a->exp+1 == a->len) {
    // Special case where no scaling is needed
    n = -a->exp;
    k = ceil((double)n*LOG) - 2;
    g = max(ceil(log(k)/log(OLDBASE)) + 1, KT);
  } else {
    // Fuck it. We need scaling by division
    fprintf(stderr, "Not implemented\n");
    abort();
  }

  bigfloat_alloc(target, k);
  convert_rec(target, k, a, n, g);

  // Carry if needed
  if(target->coef[target->len-2] > NEWBASE) {
    target->coef[target->len-1] = target->coef[target->len-2]/NEWBASE;
    target->coef[target->len-2] = target->coef[target->len-2]%NEWBASE;
  } else {
    target->len--;
  }

  // Fix exponent
  if(-a->exp+1 == a->len) {
    target->exp = -target->len+1;
  } else {
    // Fix that case later
    fprintf(stderr, "Not implemented\n");
    abort();
  }


  return;
}
