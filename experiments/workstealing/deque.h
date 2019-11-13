#ifndef DEQUE_ELEM_TYPE
typedef struct _deque_jmp_buf {
  jmp_buf self;
} deque_jmp_buf_t;

// jmp_buf: 8 registers, 16 signal masks, 1 sigmask flag
#define DEQUE_ELEM_TYPE deque_jmp_buf_t
#endif

#ifndef DEQUE_NEUTRAL_ELEM

#endif

#ifndef DEQUE_ELEM_COUNT
#define DEQUE_ELEM_COUNT (0x100)
#endif

typedef struct {
  DEQUE_ELEM_TYPE  elems[DEQUE_ELEM_COUNT]; // the elems
  size_t           top;   // points to the last filled top element
  size_t           bot;   // points to the last filled bottom element
  pthread_mutex_t  mutex; // locked when deque is accessed (implement atomically later)
} deque_t;

void deque_create(deque_t* deque);
void deque_push(deque_t* deque, DEQUE_ELEM_TYPE val);
DEQUE_ELEM_TYPE deque_pop(deque_t* deque);
void deque_unshift(deque_t* deque, DEQUE_ELEM_TYPE val);
DEQUE_ELEM_TYPE deque_shift(deque_t* deque);
