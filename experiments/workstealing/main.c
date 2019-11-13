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

/*
typedef struct _exec_ctx {
  uint64_t rip;
  uint64_t rbx;
  uint64_t rbp;
  uint64_t rdi;
  uint64_t rsi;
  uint64_t rsp;
  uint64_t r12;
  uint64_t r13;
  uint64_t r14;
  uint64_t r15;
} exec_ctx_t;
*/

#include "deque.h"

typedef struct {
  int cpu_id; // 0..N number of CPU
  deque_t work;
} thread_t;


int cpus = 0;
_Thread_local size_t cpu_id = 0;
thread_t* threads;

int spawn() {
  deque_jmp_buf_t q;
  int d = setjmp(q.self);
  printf("[%ld] setjmp returns %d\n", d);
  if(d == 0) {
    printf("[%ld] Spawning\n", cpu_id);
    deque_push(&threads[cpu_id].work, q);
    return 1;
  } else if(d == 1) {
    printf("[%ld] Returning\n", cpu_id);
  } else {
    abort();
  }

  return 0;
}

void join() {
  deque_jmp_buf_t q;
  printf("[%ld] Joining\n", cpu_id);
  q = deque_pop(&threads[cpu_id].work);
  longjmp(q.self, 1);
}

int sum(int low, int high) {
  if(high == low) {
    sleep(1); // Very expensive computation here
    return high;
  }

  int mid = low + (high-low)/2;
  printf("(%d, %d) -> (%d, %d), (%d, %d)\n", low, high, low, mid, mid+1, high);
  volatile int a, b;
  if(spawn()) {
    printf("A+\n");
    a = sum(low, mid);
  } else {
    printf("A-\n");
  }

  if(spawn()) {
    printf("B+\n");
    b = sum(mid+1, high);
  } else {
    printf("B-\n");
  }
  join();

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
  stick_core(me->cpu_id);
  cpu_id = me->cpu_id;

  printf("[%ld] CPU initialized.\n", cpu_id);

  //stall();
  sleep(1);
}

int main(void) {
  // Init N threads
  cpus = sysconf(_SC_NPROCESSORS_ONLN);
  threads  = alloca(cpus*sizeof(thread_t)); // when we return work is done and threads are dead, so this can go on the stack.
  pthread_t* pthreads = alloca(cpus*sizeof(pthread_t*));

  // Init deques
  for(size_t cpu=0;cpu<cpus;cpu++) {
    threads[cpu].cpu_id = cpu;
    deque_create(&threads[cpu].work);
  }

  // We are thread 0.
  for(size_t cpu=1;cpu<cpus;cpu++) {
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
