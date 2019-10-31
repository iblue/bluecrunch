// gcc -O2 main.c deque.c -ggdb -lpthread -o main && ./main

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h> // abort
#include <pthread.h>
#include <unistd.h> // sysconf, sleep
#include <alloca.h>
#include <malloc.h>
#include <sched.h>  // cpu_set_t, CPU_SET
#include <setjmp.h> // magic :)
#include <stdint.h>
#include <assert.h>

#include "deque.h"

int cpus = 0;

typedef struct {
  int cpu_id; // 0..N number of CPU
  void* others;
  deque_t work;
} thread_t;

int spawn() {
  // FIXME: Push execution context into deque
  return 1;
}

void stall() {
  // FIXME: Get execution context from deque and run
}

int sum(int low, int high) {
  if(high == low) {
    sleep(1); // Very expensive computation here
    return high;
  }

  int mid = low + (high-low)/2;
  printf("(%d, %d) -> (%d, %d), (%d, %d)\n", low, high, low, mid, mid+1, high);
  int a, b;
  if(spawn()) {
    a = sum(low, mid);
  }
  if(spawn()) {
    b = sum(mid+1, high);
  }
  stall();

  return a+b;
}

int stick_core(int core_id) {
  int num_cores = sysconf(_SC_NPROCESSORS_ONLN);
   if (core_id < 0 || core_id >= num_cores)
      return -1;

   cpu_set_t cpuset;
   CPU_ZERO(&cpuset);
   CPU_SET(core_id, &cpuset);

   pthread_t current_thread = pthread_self();
   return pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset);
}

void * thread_work(void *thread_ptr) {
  thread_t* me = (thread_t*)thread_ptr;

  printf("[%ld] CPU %d initialized.\n", pthread_self(), me->cpu_id);
  stick_core(me->cpu_id);

  //stall();
  sleep(1);
}

int main(void) {
  // Init N threads
  cpus = sysconf(_SC_NPROCESSORS_ONLN);
  thread_t*  threads  = alloca(cpus*sizeof(thread_t));
  pthread_t* pthreads = alloca(cpus*sizeof(pthread_t*));

  // We are thread 0.
  for(size_t cpu=1;cpu<cpus;cpu++) {
    threads[cpu].cpu_id = cpu;
    threads[cpu].others = threads;
    pthread_create(&pthreads[cpu], NULL, &thread_work, &threads[cpu]);
  }
  stick_core(0);

  // Go work.
  int total = sum(1, 100);
  printf("Done: %d\n", total);

  // Clean up
  for(size_t cpu=1;cpu<cpus;cpu++) {
    void* return_val;
    pthread_join(pthreads[cpu], &return_val);
  }
}
