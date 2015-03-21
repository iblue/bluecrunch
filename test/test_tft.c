#include "fft.h"
#include <assert.h>

complex double values[] = {6, 2, -2+2*I, 0};

int main() {
  fft_ensure_table(8); // up to 256 values
  tft_inverse(values, 3);

  assert(values[0] == 4);
  assert(values[1] == 8);
  assert(values[2] == 12);
}
