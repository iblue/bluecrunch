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

#define SIZE (1024*1024*64)
#define REP (100)

int main(void) {
  uint8_t *mem = (uint8_t*)aligned_alloc(32, SIZE);
  memset(mem, 0, SIZE);

  for(int stride=1;stride<(1024*1024);stride*=2) {
    double start = wall_clock();
    for(int i=0;i<REP;i++) {
      for(size_t i=0;i<SIZE/sizeof(uint8_t);i+=stride) {
        mem[i] += 1;
      }
    }
    double end = wall_clock();

    double time_per_loop   = (end-start)/REP;
    double bytes_processed = (double)SIZE / (double)stride;
    double gb_processed    = bytes_processed/(1024.0*1024.0*1024.0);
    double gb_per_sec      = gb_processed/time_per_loop;

    printf("Stride %d: %f GB/s [%f s]\n", stride, gb_per_sec, time_per_loop);
  }
}
