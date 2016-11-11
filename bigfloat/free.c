#include <malloc.h>
#include "bigfloat.h"

void bigfloat_free(bigfloat_t target) {
  if(target->coef != NULL) {
    free(target->coef);
    target->coef = NULL;
  }
}
