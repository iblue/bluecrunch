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

size_t int_to_str(uint32_t val, char* str) {
  for (int i=7;i>=0;i--){
    int nibble = val % 16;
    if(nibble < 10) {
      str[i] = nibble + '0';
    } else {
      str[i] = nibble - 10 + 'a';
    }

    val /= 16;
  }
  return 8;
}

size_t int_to_str10(uint32_t val, char* str) {
  for (int i=7;i>=0;i--){
    int nibble = val % 10;
    str[i] = nibble + '0';

    val /= 10;
  }
  return 8;
}

// Skipps leading zeros
size_t int_to_str_trimmed(uint32_t val, char* str) {
  if(val == 2) {
    str[0] = '2';
    return 1;
  } else {
    fprintf(stderr, "Not implemented\n");
    abort();
  }
}

size_t int_to_str_trimmed10(uint32_t val, char* str) {
  // FIXME
  if(val == 2) {
    str[0] = '2';
    return 1;
  } else {
    fprintf(stderr, "Not implemented\n");
    abort();
  }
}

// Returns length of string, fills char* string with value
size_t bigfloat_to_string(char* string, const bigfloat_t value, size_t digits, int base) {
  char* initial_string = string;
  if(value->len == 0) {
    string[0] = '0';
    string[1] = '\0';
    return 1;
  }

  int mag = value->len + value->exp;

  size_t c = value->len-1;

  if(mag == 1) {
    if(base == 10) {
      string += int_to_str_trimmed10(value->coef[c], string);
    } else {
      string += int_to_str_trimmed(value->coef[c], string);
    }
  }

  *string++ = '.';

  size_t min_c;

  if(value->len > digits/9) {
    min_c = value->len-digits/9-1;
  } else {
    min_c = 0;
  }

  if(base == 10) {
    while(c-->min_c) {
      string += int_to_str10(value->coef[c], string);
    }
  } else {
    while(c-->min_c) {
      string += int_to_str(value->coef[c], string);
    }
  }

  *string++ = '\0';

  // Truncate if required
  if(string - initial_string > digits) {
    initial_string[digits+1] = '\0';
    return digits;
  }

  return string - initial_string-1;
}

