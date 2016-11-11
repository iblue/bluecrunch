#include <malloc.h>
#include <string.h> // memcpy
#include <stdlib.h> // abort
#include <assert.h>
#include <cilk/cilk.h>    // More oomph!
#include "math.h"
#include "bigfloat.h"

#define OLDBASE (4294967296)
#define NEWBASE (100000000) // If we increase this by one digit, it will be faster (but more risk of overflow)
#define LOG (log(OLDBASE)/log(NEWBASE))
#define ALPHA (log(NEWBASE)/log(OLDBASE))
#define KT (3000) /* This seems to be optimal on my machine */

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
  double alpha = ALPHA;

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

  s->len = k;

  // Restore ptr and deallocate
  t->coef = coef_ptr;
  bigfloat_free(t);
}

#define CONVERT_DECOMPOSE_TRESH 10000
#define RADIX_CONV_TABLE_SIZE 128


size_t val_at_idx[RADIX_CONV_TABLE_SIZE] = {0};
bigfloat_t radix_conv_table[RADIX_CONV_TABLE_SIZE];

// FIXME: Clean up this mess
void ensure_radix_conversion(size_t k) {
  size_t max_k = k;
  size_t len   = 0;
  while(k >= KT) {
    k >>= 1;
    len++;
  }

  bigfloat_t current;
  bigfloat_new(current);
  bigfloat_set(current, NEWBASE);
  bigfloat_exp(current, current, k);
  bigfloat_new(radix_conv_table[0]);
  bigfloat_copy(radix_conv_table[0], current);
  printf("ensured %ld -> 0\n", k);
  val_at_idx[0] = k;

  size_t j=0;
  for(size_t i=1;i<len;i++) {
    j++;
    size_t gen_k = max_k >> (len-i);
    bigfloat_new(radix_conv_table[i]);
    bigfloat_exp(current, current, 2);
    k = gen_k;
    printf("ensured %ld -> %ld\n", k, j);
    val_at_idx[j] = k;
    bigfloat_copy(radix_conv_table[j], current);
    // FIXME: We are sometimes off by one, so we calculate that shit as well.
    // But we could do this in the radix, just taking the other value and
    // running a quick bigfloat_mului on the table entry. This eats ram.
    {
      j++;
      k += 1;
      bigfloat_realloc(current, current->len+1); // Prevent overflow
      bigfloat_mului(current, NEWBASE);
      printf("ensured %ld -> %ld\n", k, j);
      val_at_idx[j] = k;
      bigfloat_copy(radix_conv_table[j], current);
    }

    // FIXME: Cache the FFT transforms where needed
  }

  // Terminator
  bigfloat_new(radix_conv_table[len]);
  radix_conv_table[len]->coef = NULL;

  bigfloat_free(current);
}

size_t exp_to_radix_tbl_entry(size_t k) {
  for(size_t i=0;i<RADIX_CONV_TABLE_SIZE;i++) {
    if(val_at_idx[i] == k) {
      return i;
    }
  }

  fprintf(stderr, "did not find %ld in tbl\n", k);
  abort();
}

void free_radix_conversion() {
  size_t i = 0;
  while(1) {
    bigfloat_free(radix_conv_table[i]);
    i++;
    if(radix_conv_table[i]->coef == NULL) {
      break;
    }
  }
}

void convert_rec(bigfloat_t s, size_t k, const bigfloat_t y, size_t n, size_t g) {
  if(k <= KT) {
    bigfloat_alloc(s, k);
    convert_trunc(s, y, k, n);
  } else {
    size_t kh = (k+1)/2;
    size_t kl = k - kh + 1;

    // Choose nh such that 4gb^kh < 2^nh
    size_t nh = floor((double)(kh*log((double)NEWBASE)+log(4.)+log(g))/log((double)OLDBASE)+1);

    // Choose nl such that 4gb^kl < 2^nl
    size_t nl = floor((double)(kl*log((double)NEWBASE)+log(4.)+log(g))/log((double)OLDBASE)+1);

    // yh = floor(y*(2^32)^(nh-n)) = y/[(2^32)^(n-nh)]
    bigfloat_t yh;
    bigfloat_new(yh);
    bigfloat_alloc(yh, y->len-(n-nh));
    memcpy(yh->coef, y->coef+(n-nh), yh->len*sizeof(uint32_t));

    // yl = b^(k-kl)*y mod (2^32)^n bdiv (2^32)^(n-nl)
    bigfloat_t yl;
    bigfloat_new(yl);
    // TODO: Performance: This can be done as middle product
    printf("using table for %ld\n", k-kl);
    bigfloat_mul(yl, y, radix_conv_table[exp_to_radix_tbl_entry(k-kl)], 0);
    yl->len = n; // yl <- yl mod (2^32)^n
    uint32_t *yl_coef_ptr = yl->coef;
    yl->coef += (n - nl); // bdiv (2^32)^(n-nl)
    yl->len  -= (n - nl); // vvvvvvvvvvvvvvvvv

    // Recurse
    bigfloat_t sh;
    bigfloat_new(sh);
    bigfloat_t sl;
    bigfloat_new(sl);

    // Run parallel
    if(k > CONVERT_DECOMPOSE_TRESH) {
      cilk_spawn convert_rec(sh, kh, yh, nh, g);
      convert_rec(sl, kl, yl, nl, g);
      cilk_sync;
    } else {
      convert_rec(sh, kh, yh, nh, g);
      convert_rec(sl, kl, yl, nl, g);
    }

    // fixups. if the trailing digit of sh is b-1 and the leading digit of sl is 0
    if(sh->coef[0] == NEWBASE-1 && sl->coef[sl->len-1] == 0) {
      // add 1. FIXME: Move into a bigfloat_addu()
      uint32_t carry = 1;
      for(size_t i=0;i<sh->len;i++) {
        uint64_t sum = sh->coef[i] + carry;
        uint32_t upper = sum / NEWBASE;
        uint32_t lower = sum % NEWBASE;
        carry = upper;
        sh->coef[i] = lower;

        if(carry == 0) {
          break;
        }
      }
    }

    // fixup. if the trailing digit of sh is 0 and the leading digit of sl is b-1
    char sl_is_zero = 0;
    if(sh->coef[0] == 0 && sl->coef[sl->len-1] == NEWBASE-1) {
      // sl = 0
      sl_is_zero = 1;
    }

    // s = floor(sh/b)*b^(kl) + sl (s is now in base b)
    bigfloat_alloc(s, sh->len-1+kl);
    memset(s->coef, 0, kl*sizeof(uint32_t)); // fill unitialized mem
    memcpy(s->coef+kl, sh->coef+1, (sh->len-1)*sizeof(uint32_t)); // skip the first digit (this is equiv to div by 2^32 and flooring)
    if(!sl_is_zero) {
      // FIXME: Do we really need to add? There should be no carrying, so we
      // can just copy. Test this.
      bigfloat_add(s, s, sl, 0);
    }

    bigfloat_free(sh);
    bigfloat_free(sl);

    bigfloat_free(yh);
    yl->coef = yl_coef_ptr;
    bigfloat_free(yl);
  }
}

// Converts from OLDBASE to NEWBASE
void bigfloat_radix(bigfloat_t target, const bigfloat_t a) {
  size_t n, k, g;
  if(-a->exp+1 == a->len) {
    // Special case where no scaling is needed
    n = -a->exp;
    k = ceil((double)n*LOG) - 1;
    g = max(ceil(log(k)/log(OLDBASE)) + 1, KT);
  } else {
    // Fuck it. We need scaling by division
    fprintf(stderr, "Not implemented\n");
    abort();
  }

  ensure_radix_conversion(k);
  bigfloat_alloc(target, k);
  convert_rec(target, k, a, n, g);
  free_radix_conversion();

  // Carry if needed
  if(target->coef[target->len-1] > NEWBASE) {
    bigfloat_realloc(target, target->len+1);
    target->coef[target->len-1] = target->coef[target->len-2]/NEWBASE;
    target->coef[target->len-2] = target->coef[target->len-2]%NEWBASE;
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
