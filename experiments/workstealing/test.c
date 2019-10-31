// gcc test.c deque.c -ggdb -O2 -DDEQUE_ELEM_COUNT=5 -DDEQUE_ELEM_TYPE=uint32_t -lpthread -o test && ./test

#include <stdio.h>
#include <stdlib.h> // abort
#include <pthread.h>
#include <stdint.h>
#include <assert.h>

#include "deque.h"

int main(void) {
  // Test deque.
  deque_t magic[1];
  deque_create(magic);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1); // [1]
  assert(deque_pop(magic) == 1);  // []
  assert(deque_pop(magic) == 0); // []
  deque_push(magic, 1); // [1]
  deque_push(magic, 2); // [1, 2]
  deque_push(magic, 3); // [1, 2, 3]
  deque_push(magic, 4); // [1, 2, 3, 4]
  assert(deque_pop(magic) == 4); // [1, 2, 3]
  assert(deque_pop(magic) == 3); // [1, 2]
  assert(deque_pop(magic) == 2); // [1]
  assert(deque_pop(magic) == 1); // []
  assert(deque_pop(magic) == 0); // []
  deque_push(magic, 1002);  // [1002]
  deque_push(magic, 1003);  // [1002, 1003]
  deque_unshift(magic, 1001); // [1001, 1002, 1003]
  assert(deque_pop(magic) == 1003); // [1001, 1002]
  assert(deque_pop(magic) == 1002); // [1001]
  assert(deque_pop(magic) == 1001); // []
  assert(deque_pop(magic) == 0);
  deque_push(magic, 1202); // [1202]
  deque_push(magic, 1203); // [1202, 1203]
  deque_unshift(magic, 1201); // [1201, 1202, 1203]
  assert(deque_pop(magic) == 1203); // ...
  assert(deque_pop(magic) == 1202);
  assert(deque_pop(magic) == 1201);
  assert(deque_pop(magic) == 0);
  deque_push(magic, 2002);  // [2002]
  deque_push(magic, 2003); // [2002, 2003]
  deque_unshift(magic, 2001); // [2001, 2002, 2003]
  assert(deque_shift(magic) == 2001); // [2002, 2003]
  assert(deque_shift(magic) == 2002); // [2003]
  assert(deque_shift(magic) == 2003); // []
  deque_unshift(magic, 2004); // [2004]
  assert(deque_shift(magic) == 2004); // []
  deque_push(magic, 100); // [100]
  deque_push(magic, 101); // [100, 101]
  deque_push(magic, 102); // [100, 101, 102]
  deque_push(magic, 103); // [100, 101, 102, 103]
  assert(deque_shift(magic) == 100); // [101, 102, 103]
  deque_push(magic, 104); // [101, 102, 103, 104]
  assert(deque_shift(magic) == 101); // [102, 103, 104]
  assert(deque_pop(magic) == 104); // [102, 103]
  deque_push(magic, 105); // [102, 103, 105]
  assert(deque_shift(magic) == 102); // [103, 105]
  assert(deque_pop(magic) == 105); // [103]
  assert(deque_shift(magic) == 103); // []
  assert(deque_pop(magic) == 0);

  return 0;
}
