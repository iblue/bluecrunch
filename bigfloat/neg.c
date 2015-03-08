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

void bigfloat_neg(bigfloat_t num) {
  if(bigfloat_iszero(num)) {
    return;
  }

  num->sign = !num->sign;
}
