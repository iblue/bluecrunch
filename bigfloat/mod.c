//  SIMD
#include <malloc.h>
#include <pmmintrin.h>
#include <string.h>

#include <cilk/cilk.h>
#include "bigfloat.h"
#include "fft.h"

// Takes modulo a 64 bit unsigned integer
uint64_t bigfloat_modu(bigfloat_t a, uint64_t m) {
  uint64_t base = 1;
  uint64_t result = 0;
  for(size_t i=0;i<a->len;i++) {
    uint64_t v = a->coef[i] % m;
    __uint128_t product = (__uint128_t)v*(__uint128_t)base;
    product %= m;
    __uint128_t newbase = (__uint128_t)base * (__uint128_t)m;
    newbase  %= m;
    __uint128_t sum = (__uint128_t)result + product;
    result += sum;
    base    = newbase;
  }

  return result;
}
