#include <stdio.h>
#include <assert.h>
#include <complex.h>
#include "fft.h"

int main(void) {
  // 12 bits, 1/3 words
  {
    complex double T[16] = {0xbcd, 0x34a, 0x12};

    // Scale
    for(size_t i=0;i<16;i++) {
      T[i] *= 16;
    }

    uint32_t words[1];

    fft_to_int(T, 4, words, 1, 12);

    assert(words[0] == 0x1234abcd);
  }

  // 12 bits, 2/3 words
  {
    complex double T[16] = {0xbcd, 0x34a, 0xa12, 0xbba, 0xdcc, 0xd};

    // Scale
    for(size_t i=0;i<16;i++) {
      T[i] *= 16;
    }

    uint32_t words[2];

    fft_to_int(T, 4, words, 2, 12);

    assert(words[0] == 0x1234abcd);
    assert(words[1] == 0xddccbbaa);
  }

  // 12 bits, 2/3 words
  {
    complex double T[16] = {0xbcd, 0x34a, 0xa12, 0xbba, 0xdcc, 0x11d, 0x322, 0x443};

    // Scale
    for(size_t i=0;i<16;i++) {
      T[i] *= 16;
    }

    uint32_t words[3];

    fft_to_int(T, 4, words, 3, 12);

    assert(words[0] == 0x1234abcd);
    assert(words[1] == 0xddccbbaa);
    assert(words[2] == 0x44332211);
  }

  return 0;
}
