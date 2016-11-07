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

// Calculates approximate logarithm
double bigfloat_log(bigfloat_t target) {
  // t =~ target->coef[target->len-1] * (2^32)^(target->exp + target->len-1) +
  //      target->coef[target->len-2] * (2^32)^(target->exp + target->len-2) // if len > 1
  //
  double approx = 0.0;

  if(target->len > 0) {
    approx = target->coef[target->len-1] * pow(4294967296.0, target->exp + target->len-1);
  }
  if(target->len > 1) {
    approx += target->coef[target->len-1] * pow(4294967296.0, target->exp + target->len-2);
  }

  return log(approx);
}
