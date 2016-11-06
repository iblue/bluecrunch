#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>

//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <omp.h>
#include "fft.h"
#include "bigfloat.h"

void bigfloat_exp(bigfloat_t target, const bigfloat_t a, uint64_t exp, int tds) {
  bigfloat_t cp; // current potency
  bigfloat_new(cp);
  bigfloat_copy(cp, a);

  bigfloat_set(target, 1);

  while(exp > 0) {
    if(exp % 2 == 1) {
      bigfloat_mul(target, target, cp, 0, tds);
    }

    bigfloat_mul(cp, cp, cp, 0, tds);

    exp /= 2;
  }
}
