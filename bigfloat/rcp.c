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

void bigfloat_rcp(bigfloat_t target, const bigfloat_t a, size_t p, int tds) {
  //  Compute reciprocal using Newton's Method.

  //  r1 = r0 - (r0 * x - 1) * r0

  if (a->len == 0) {
    fprintf(stderr, "Divide by Zero\n");
    abort();
  }

  //  Collect operand
  int64_t Aexp = a->exp;
  size_t AL = a->len;
  uint32_t *AT = a->coef;

  //  End of recursion. Generate starting point.
  if (p == 0){
      //  Truncate precision to 3.
      p = 3;
      if (AL > p){
          size_t chop = AL - p;
          AL = p;
          Aexp += chop;
          AT += chop;
      }

      //  Convert number to floating-point.
      double val = AT[0];
      if (AL >= 2)
          val += AT[1] * 1000000000.;
      if (AL >= 3)
          val += AT[2] * 1000000000000000000.;

      //  Compute reciprocal.
      val = 1. / val;
      Aexp = -Aexp;

      //  Scale
      while (val < 1000000000.){
          val *= 1000000000.;
          Aexp--;
      }

      //  Rebuild a bigfloat_t.
      uint64_t val64 = (uint64_t)val;

      target->sign = a->sign;

      target->coef = (uint32_t*)malloc(sizeof(uint32_t)*2);
      target->coef[0] = (uint32_t)(val64 % 1000000000);
      target->coef[1] = (uint32_t)(val64 / 1000000000);
      target->len = 2;
      target->exp = Aexp;

      return;
  }

  //  Half the precision
  size_t s = p / 2 + 1;
  if (p == 1) s = 0;
  if (p == 2) s = 1;

  //  Recurse at half the precision
  bigfloat_t T;
  bigfloat_new(T);
  bigfloat_rcp(T, a, s, tds);

  //  r1 = r0 - (r0 * x - 1) * r0
  bigfloat_t one;
  bigfloat_new(one);
  bigfloat_set(one, 1);
  bigfloat_t tmp;
  bigfloat_new(tmp);
  bigfloat_mul(tmp, a, T, p, tds);
  bigfloat_t tmp2;
  bigfloat_new(tmp2);
  bigfloat_sub(tmp2, tmp, one, p);
  bigfloat_t tmp3;
  bigfloat_new(tmp3);
  bigfloat_mul(tmp3, tmp2, T, p, tds);
  bigfloat_free(tmp); // because overwrite
  bigfloat_sub(target, T, tmp3, p);
  bigfloat_free(tmp3);
  bigfloat_free(T);
  bigfloat_free(one);
}

