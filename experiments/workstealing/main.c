#include <stdio.h>
#include <pthread.h>
#include <unistd.h> // sysconf

int sum(int low, int high) {

  if(high == low) {
    return high;
  }

  int mid = low + (high-low)/2;
  printf("(%d, %d) -> (%d, %d), (%d, %d)\n", low, high, low, mid, mid+1, high);
  int a = sum(low, mid);
  int b = sum(mid+1, high);

  return a+b;
}


int cpus=0;

int main(void) {
  cpus = sysconf(_SC_NPROCESSORS_ONLN);
  printf("cores: %d\n", cpus);
  printf("%d\n", sum(0, 100));
}
