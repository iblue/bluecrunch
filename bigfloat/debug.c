#include <stdint.h>
#include <malloc.h>
#include "bigfloat.h"

void bigfloat_print(const char* name, const bigfloat_t value) {
  printf("%s = bigfloat_t {exp = %ld, len = %ld, sign = %d, coef = ",
      name, value->exp, value->len, value->sign);
  if(value->coef == NULL) {
    printf("NULL");
  } else {
    printf("{ ");
    for(size_t i=0; i<value->len; i++) {
      printf("%d", value->coef[i]);
      if(i != value->len - 1) {
        printf(", ");
      }
    }
    printf(" }");
  }
  printf("}\n");
}
