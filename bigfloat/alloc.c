#include <malloc.h>
#include <string.h> // memset
#include "bigfloat.h"

void bigfloat_alloc(bigfloat_t target, size_t size) {
  bigfloat_free(target);

  target->exp     = 0;
  target->sign    = 1;
  target->len     = size;
  target->coef    = (uint32_t*) malloc(target->len*sizeof(uint32_t));
}

void bigfloat_realloc(bigfloat_t target, size_t size) {
  size_t old_size = target->len;
  target->len     = size;
  target->coef    = (uint32_t*) realloc(target->coef, target->len*sizeof(uint32_t));

  // Zero new data.
  if(target->len > old_size) {
    memset(target->coef + old_size, 0, (target->len - old_size)*sizeof(uint32_t));
  }
}
