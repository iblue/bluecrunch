#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 700
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <time.h>
#include "bench.h"

inline double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

double time_start = 0;

void task_start(char *description) {
  printf("%s...", description);
  fflush(stdout);
  time_start = wall_clock();
}

void task_end() {
  double time_end = wall_clock();

  printf(" ok [%f seconds]\n", time_end - time_start);
}
