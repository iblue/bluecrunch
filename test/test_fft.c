#include "fft.h"
#include <assert.h>


int main() {
  fft_ensure_table(8); // up to 256 values

  {
    __attribute__ ((aligned (32))) complex double values[] = {1, 2, 3, 4};

    fft_forward(values, 2, 1);

    // FIXME: FP accuracy!
    /*
    assert(values[0] == 10);
    assert(values[1] == -2);
    assert(values[2] == -2-2*I);
    assert(values[3] == -2+2*I);
    */
  }
}
