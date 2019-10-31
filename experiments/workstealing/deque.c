#include "deque.h"

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


typedef struct {
  DEQUE_ELEM_TYPE* elems; // the elems
  size_t           top;   // points to the last filled top element
  size_t           bot;   // points to the last filled bottom element
  pthread_mutex_t* mutex; // locked when deque is accessed (implement atomically later)
} deque_t;

void deque_create(deque_t* deque) {
  deque->mutex = malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(deque->mutex, NULL);
  pthread_mutex_lock(deque->mutex);
  deque->elems = malloc(DEQUE_ELEM_COUNT*sizeof(uint32_t));
  deque->size = size;
  deque->top = 0;
  deque->bot = 0;
  pthread_mutex_unlock(deque->mutex);
}

void deque_push(deque_t* deque, DEQUE_ELEM_TYPE val) {
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

DEQUE_ELEM_TYPE deque_pop(deque_t* deque) {
  pthread_mutex_lock(deque->mutex);
  DEQUE_ELEM_TYPE x;
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

void deque_unshift(deque_t* deque, DEQUE_ELEM_TYPE val) {
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

DEQUE_ELEM_TYPE deque_shift(deque_t* deque) {
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
