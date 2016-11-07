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

void bigfloat_div(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds) {
  //  Division
  bigfloat_t rcp;
  bigfloat_new(rcp);

  if(p == 0) {
    //  Default value
    p = a->len + b->len;
  }

  bigfloat_rcp(rcp, b, p, tds);
  bigfloat_mul(target, a, rcp, p, tds);
  bigfloat_free(rcp);
}
