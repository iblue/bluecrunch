#define DEQUE_ELEM_TYPE uint32_t
#define DEQUE_ELEM_SIZE (5)

#include "deque.h"

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
}
