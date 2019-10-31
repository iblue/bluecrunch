// gcc -O2 main.c -ggdb -lpthread -o main && ./main

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

int cpus = 0;

typedef struct {
  uint32_t* elems;        // the elems
  uint32_t  size;         // deque size in elems
  uint32_t  top;          // points to the last filled top element
  uint32_t  bot;          // points to the last filled bottom element
  pthread_mutex_t* mutex; // locked when deque is accessed (implement atomically later)
} deque_t;

void deque_create(deque_t* deque, uint32_t size) {
  deque->mutex = malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(deque->mutex, NULL);
  pthread_mutex_lock(deque->mutex);
  deque->elems = malloc(size*sizeof(uint32_t));
  deque->size = size;
  deque->top = 0;
  deque->bot = 0;
  pthread_mutex_unlock(deque->mutex);
}

void deque_push(deque_t* deque, uint32_t val) {
  pthread_mutex_lock(deque->mutex);
  deque->top++;
  if(deque->top >= deque->size) { // overflow
    deque->top = 0;
  }
  if(deque->top == deque->bot) {
    printf("deque_push: Deque Full\n");
    abort();
  }
  deque->elems[deque->top] = val;
  pthread_mutex_unlock(deque->mutex);
}

uint32_t deque_pop(deque_t* deque) {
  pthread_mutex_lock(deque->mutex);
  uint32_t x;
  if(deque->top != deque->bot) {
    x = deque->elems[deque->top];
    deque->top--;
    if(deque->top >= deque->size) { // underflow
      deque->top = deque->size-1;
    }
  } else {
    x = 0;
  }
  pthread_mutex_unlock(deque->mutex);
  return x;
}

void deque_unshift(deque_t* deque, uint32_t val) {
  pthread_mutex_lock(deque->mutex);
  deque->elems[deque->bot] = val;
  deque->bot--;
  if(deque->bot >= deque->size) { // underflow
    deque->bot = deque->size-1;
  }
  if(deque->top == deque->bot) {
    printf("deque_unshift: Deque Full\n");
    abort();
  }
  pthread_mutex_unlock(deque->mutex);
}

uint32_t deque_shift(deque_t* deque) {
  pthread_mutex_lock(deque->mutex);
  uint32_t x;
  if(deque->top != deque->bot) {
    deque->bot++;
    if(deque->bot >= deque->size) { // overflow
      deque->bot = 0;
    }
    x = deque->elems[deque->bot];
  } else {
    x = 0;
  }
  pthread_mutex_unlock(deque->mutex);
  return x;
}

typedef struct {
  int cpu_id; // 0..N number of CPU
  void* others;
  deque_t work;
} thread_t;

int spawn() {
  // if current thread, return 1
  //return !setjmp(buf);
}

void stall() {

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
  // Test deque.
  deque_t magic[1];
  deque_create(magic, 5);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1);
  assert(deque_pop(magic) == 1);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1);
  deque_push(magic, 2);
  deque_push(magic, 3);
  deque_push(magic, 4);
  assert(deque_pop(magic) == 4);
  assert(deque_pop(magic) == 3);
  assert(deque_pop(magic) == 2);
  assert(deque_pop(magic) == 1);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1002);
  deque_push(magic, 1003);
  deque_unshift(magic, 1001);
  assert(deque_pop(magic) == 1003);
  assert(deque_pop(magic) == 1002);
  assert(deque_pop(magic) == 1001);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1202);
  deque_push(magic, 1203);
  deque_unshift(magic, 1201);
  assert(deque_pop(magic) == 1203);
  assert(deque_pop(magic) == 1202);
  assert(deque_pop(magic) == 1201);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 2002);
  deque_push(magic, 2003);
  deque_unshift(magic, 2001);
  assert(deque_shift(magic) == 2001);
  assert(deque_shift(magic) == 2002);
  assert(deque_shift(magic) == 2003);
  deque_unshift(magic, 2004);
  assert(deque_shift(magic) == 2004);
  deque_push(magic, 100);
  deque_push(magic, 100);
  deque_push(magic, 100);
  deque_push(magic, 100);
  deque_push(magic, 100);
  assert(deque_shift(magic) == 0);
  return 1;


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
