#include <stdio.h>
#include <pthread.h>
#include <cilk/cilk.h>
#include <unistd.h> // sysconf, sleep

int sum(int low, int high) {
  if(high == low) {
    sleep(1); // Very expensive computation here
    return high;
  }

  int mid = low + (high-low)/2;
  printf("(%d, %d) -> (%d, %d), (%d, %d)\n", low, high, low, mid, mid+1, high);
  int a = cilk_spawn sum(low, mid);
  int b = sum(mid+1, high);

  cilk_sync;

  return a+b;
}


int cpus=0;

int main(void) {
  cpus = sysconf(_SC_NPROCESSORS_ONLN);
  printf("cores: %d\n", cpus);
  printf("%d\n", sum(0, 100));
}
