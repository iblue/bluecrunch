#include <stdio.h>

__attribute__ ((noinline)) size_t _table_select(size_t length) {
  size_t p     = __builtin_popcountl(length); // 1 for 2*2^k, 2 for 3*2^k
  size_t k     = __builtin_ctzl(length); // returns k or k+1 for 2*2^k
  p--; // 0 for 2*2^k, 1 for 3*2^k
  k -= p ? 0 : 1; // correct k by -1 if 2*2^k
  size_t base  = 1 << k; // 2^k
  size_t ret   = base << 1; // 2*2^k
  ret += p ? base : 0; // adds 2^k if p == 1 (leading 2*2^k if p == 0, 3*2^k else)

  // did we get it right? => return
  if(ret == length) {
    return (k<<1)+p;
  } else {
    return -1;
  }
}

int main(void) {
  printf("%ld\n", _table_select(64));
  return 0;
}
