#include <stdlib.h> // abort
#include <stdio.h> // printf
#include <pthread.h>
#include <stdint.h>
#include <assert.h>
#include <setjmp.h> // jmp_buf

#include "deque.h"

/**
 * Double ended queue
 *
 * * Single producer
 * * Many consumers
 * * Lock free (todo)
 * * size fixed at compile time
 * * element size > 128 bits
 * * thread- and interrupt safe
 */

void deque_create(deque_t* deque) {
  pthread_mutex_init(&deque->mutex, NULL);
  pthread_mutex_lock(&deque->mutex);
  deque->top = 0;
  deque->bot = 0;
  pthread_mutex_unlock(&deque->mutex);
}

void deque_push(deque_t* deque, DEQUE_ELEM_TYPE val) {
  pthread_mutex_lock(&deque->mutex);
  deque->top++;
  if(deque->top >= DEQUE_ELEM_COUNT) { // overflow
    deque->top = 0;
  }
  if(deque->top == deque->bot) {
    printf("deque_push: Deque Full\n");
    abort();
  }
  deque->elems[deque->top] = val;
  pthread_mutex_unlock(&deque->mutex);
}

DEQUE_ELEM_TYPE deque_pop(deque_t* deque) {
  pthread_mutex_lock(&deque->mutex);
  DEQUE_ELEM_TYPE x;
  if(deque->top != deque->bot) {
    x = deque->elems[deque->top];
    deque->top--;
    if(deque->top >= DEQUE_ELEM_COUNT) { // underflow
      deque->top = DEQUE_ELEM_COUNT-1;
    }
  } else {
    // FIXME: Return neutral element
    // x = NULL;
    printf("Queue empty\n");
    abort();
  }
  pthread_mutex_unlock(&deque->mutex);
  return x;
}

void deque_unshift(deque_t* deque, DEQUE_ELEM_TYPE val) {
  pthread_mutex_lock(&deque->mutex);
  deque->elems[deque->bot] = val;
  deque->bot--;
  if(deque->bot >= DEQUE_ELEM_COUNT) { // underflow
    deque->bot = DEQUE_ELEM_COUNT-1;
  }
  if(deque->top == deque->bot) {
    printf("deque_unshift: Deque Full\n");
    abort();
  }
  pthread_mutex_unlock(&deque->mutex);
}

DEQUE_ELEM_TYPE deque_shift(deque_t* deque) {
  pthread_mutex_lock(&deque->mutex);
  DEQUE_ELEM_TYPE x;
  if(deque->top != deque->bot) {
    deque->bot++;
    if(deque->bot >= DEQUE_ELEM_COUNT) { // overflow
      deque->bot = 0;
    }
    x = deque->elems[deque->bot];
  } else {
    // FIXME: Return neutral element
    // x = NULL;
    printf("Queue empty\n");
    abort();
  }
  pthread_mutex_unlock(&deque->mutex);
  return x;
}
