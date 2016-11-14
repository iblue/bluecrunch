#if __STDC_VERSION__ >= 199901L
#define _XOPEN_SOURCE 700
#else
#define _XOPEN_SOURCE 500
#endif /* __STDC_VERSION__ */

#include <stdio.h>
#include <time.h>
#include "bench.h"

double wall_clock() {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1.0e-9*t.tv_nsec;
}

#define MAX_LEVELS 3

int last_level = 0;
double time_start[MAX_LEVELS];

void task_start(int level, char *description) {
  if(last_level == level-1) {
    printf("\n");
  }
  printf("%*s%s...", 2*level, "", description);
  fflush(stdout);
  time_start[level] = wall_clock();
  last_level = level;
}

void task_end(int level) {
  double time_end = wall_clock();

  if(last_level != level) {
    printf("%*sDone [%f seconds]\n", 2*level, "", time_end - time_start[level]);
  } else {
    printf(" done [%f seconds]\n", time_end - time_start[level]);
  }
  time_start[level] = 0;
}
