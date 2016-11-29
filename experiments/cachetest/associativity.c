#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 700
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

#define SIZE (1024*16)
#define REP (1)

// Bulldozer:
// L1: 16K 4-way associative with 64 Byte Cache lines (per Core) (256 lines)
// L2:  2M 16-way associative with 64 Byte Cache lines (shared among 2 Cores) (32768 lines)
// L3:  8M 64-way associative with 64 Byte Cache lines (shared among all Cores) (131072 lines)
// 1 Cache Line = 64 Bytes = 8 doubles = 4 complex doubles
// To use all lines, we need 64 B Padding in each row.
// AVX registers: 512 bytes = 8 Cache Lines

int main(void) {
  uint8_t *mem = (uint8_t*)aligned_alloc(32, SIZE);
  memset(mem, 0, SIZE);

  for(int stride=1;stride<SIZE;stride*=2) {
    {
      double start = wall_clock();
      for(int rep=0;rep<REP;rep++) {
        for(int j=0;j<stride;j++) {
          for(int i=0;i<SIZE/stride;i++) {
            //printf("(%d,%d) -> %d\n", i, j, i*stride+j);
            mem[i*stride+j] += 1;
          }
        }
      }
      double end = wall_clock();

      double time_per_loop   = (end-start)/REP;
      printf("Stride %d: %f s (unpadded)\n", stride, time_per_loop);
    }
  }
}
