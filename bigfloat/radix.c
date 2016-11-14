#include <malloc.h>
#include <string.h> // memcpy
#include <stdlib.h> // abort
#include <assert.h>
#include <cilk/cilk.h>    // More oomph!
#include "bench.h"
#include "math.h"
#include "bigfloat.h"

#define OLDBASE (4294967296)
#define NEWBASE (1000000000) // If we increase this by one digit, it will be faster (but more risk of overflow)
#define LOG (log(OLDBASE)/log(NEWBASE))
#define ALPHA (log(NEWBASE)/log(OLDBASE))
#define KT (1000) /* This seems to be optimal on my machine */

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
bigfloat_t radix_conv_table[128];

size_t find(size_t k) {
  for(size_t i=0;i<RADIX_CONV_TABLE_SIZE;i++) {
    if(val_at_idx[i] == k) {
      return i;
    }
  }

  return -1;
}

void insert_unless_exists(size_t k) {
  for(size_t i=0;i<RADIX_CONV_TABLE_SIZE;i++) {
    if(val_at_idx[i] == k) {
      break;
    }
    if(val_at_idx[i] == 0) {
      val_at_idx[i] = k;
      break;
    }
  }
}

int compare(const void *a, const void* b) {
  return (int)(*(size_t*)a - *(size_t*)b);
}

void sort() {
  size_t len = 0;

  // Calculate length
  for(size_t i=0;i<RADIX_CONV_TABLE_SIZE;i++) {
    if(val_at_idx[i] == 0) {
      len = i-1;
      break;
    }
  }

  qsort(val_at_idx, len, sizeof(size_t), compare);
}

void build_table(size_t k) {
  if(k < KT) {
    return;
  }

  size_t kh  = (k+1)/2;
  size_t kl  = k - kh + 1;
  size_t exp = k - kl;

  build_table(kh);
  if(kl != kh) {
    build_table(kl);
  }

  insert_unless_exists(exp);
}

void generate_factors() {
  // Generate first
  bigfloat_new(radix_conv_table[0]);
  bigfloat_set(radix_conv_table[0], NEWBASE);
  bigfloat_exp(radix_conv_table[0], radix_conv_table[0], val_at_idx[0]);

  //printf("ensured %ld -> 0 by inital\n", val_at_idx[0]);

  for(size_t i=1;i<RADIX_CONV_TABLE_SIZE;i++) {
    //printf("processing %ld -> %ld\n", val_at_idx[i], i);
    if(val_at_idx[i] == 0) {
      break;
    }

    size_t current_exp = val_at_idx[i];

    // Check if we find exp - 1
    size_t idx = find(current_exp - 1);
    if(idx != -1) {
      // Simple job. We have B^(exp-1). Just multiply by B.
      bigfloat_new(radix_conv_table[i]);
      bigfloat_alloc(radix_conv_table[i], radix_conv_table[idx]->len+1);
      bigfloat_mulu(radix_conv_table[i], radix_conv_table[idx], NEWBASE);
      //printf("ensured %ld -> %ld by multiplication\n", current_exp, i);
      continue;
    }

    // Check if we find exp/2
    idx = find(current_exp/2);
    if(idx != -1) {
      // Ok, just square
      bigfloat_new(radix_conv_table[i]);
      bigfloat_mul(radix_conv_table[i], radix_conv_table[idx], radix_conv_table[idx], 0);

      // New, if exp is odd, we are missing a multiplication. We have
      // 2*(exp/2), we need 2*(exp/2)+2.
      if((current_exp & 1) == 1) {
        bigfloat_realloc(radix_conv_table[i], radix_conv_table[i]->len+1);
        bigfloat_mului(radix_conv_table[i], NEWBASE);
        //printf("ensured %ld -> %ld by squaring and multiplication\n", current_exp, i);
      } else {
        //printf("ensured %ld -> %ld by squaring\n", current_exp, i);
      }
      continue; // Either way, we're done
    }

    fprintf(stderr, "did not find previous value for %ld in tbl\n", current_exp);
    abort();
  }
}

void ensure_radix_conversion(size_t k) {
  if(k < KT) {
    return;
  }

  //printf("Ensuring table (k = %ld)\n", k);
  build_table(k);
  sort();
  generate_factors();
}

void free_radix_conversion() {
  for(size_t i=0;i<RADIX_CONV_TABLE_SIZE;i++) {
    bigfloat_free(radix_conv_table[i]);

    if(val_at_idx[i] == 0) {
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
    //printf("using table for %ld\n", k-kl);
    bigfloat_mul(yl, y, radix_conv_table[find(k-kl)], 0);

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

  task_start(2, "Precomputing base conversion table");
  ensure_radix_conversion(k);
  task_end(2);
  bigfloat_alloc(target, k);
  task_start(2, "Conversion");
  convert_rec(target, k, a, n, g);
  task_end(2);
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
