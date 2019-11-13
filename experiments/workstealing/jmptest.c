#include <setjmp.h>
#include <stdio.h>

jmp_buf me;

void dead() {
}

int main(void) {
  volatile int a = 0;
  int b = 0;
  if(!setjmp(me)) {
    a = 1;
    b = 1;
    printf("first\n");
    printf("a = %d\n", a);
    printf("b = %d\n", b);
    a = 2;
    b = 2;
    longjmp(me, 1);
    goto dead;
  } else {
    printf("second\n");
    printf("a = %d\n", a);
    printf("b = %d\n", b);
  }

  return 0;
dead:
    printf("a = %d\n", a);
    printf("b = %d\n", b);
    return 0;
}
