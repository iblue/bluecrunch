/*
gcc -std=c11 -Wall -Werror -ggdb -O3 -mavx -fcilkplus ./test.c -o magic
rm callgrind.out.*
valgrind --tool=callgrind --dump-instr=yes ./magic
kcachegrind callgrind.out.*
*/

#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>

#define CHECKPRIME (0x9f063b05)

uint32_t checksum(uint32_t *AT, size_t AL) {
  uint64_t base = 1;
  uint64_t prime = CHECKPRIME;
  uint64_t sum = 0;
  for(size_t i=0;i<AL;i++) {
    uint64_t v = AT[i];
    uint64_t product = (v*base) % prime;
    base = (base << 32)%prime;
    sum = sum + product;
    sum %= prime;
  }

  return (uint32_t)sum;
}

uint32_t checksum_mul(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p  = (av*bv)%CHECKPRIME;
  return (uint32_t)p;
}

uint32_t checksum_add(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p  = (av+bv)%CHECKPRIME;
  return (uint32_t)p;
}

uint32_t checksum_sub(uint32_t a, uint32_t b) {
  uint64_t av = a;
  uint64_t bv = b;
  uint64_t p;
  if(av > bv) {
    p = av-bv;
  } else {
    p = av+CHECKPRIME-bv;
  }
  p  = p%CHECKPRIME;
  return (uint32_t)p;
}


__attribute__ ((noinline)) void _ll_sub_inplace(uint32_t *CT, size_t CL, uint32_t *AT, size_t AL) {
  uint64_t* CT64 = (uint64_t*) CT;
  size_t    CL64 = CL/2;
  uint64_t* AT64 = (uint64_t*) AT;
  size_t    AL64 = AL/2;

  int carry = 0;
  size_t i64;
  size_t i32;

  // Operate on 64 bits first
  for(i64 = 0;i64<AL64;i64++) {
    uint64_t a = CT64[i64];
    uint64_t b = AT64[i64];
    uint64_t c = a - b - carry;
    carry = (c>a);
    CT64[i64] = c;
  }

  // Now, if there is an odd number of coefs in AL, we need to handle one 32
  // bit sub
  if(AL%2 == 1) {
    i32 = i64*2;

    uint32_t a = CT[i32];
    uint32_t b = AT[i32];
    uint32_t c = a - b - carry;
    carry = (c>a);
    CT[i32] = c;
    i32++;

    // Now, if there is more in CT we need to handle another 32 bit addition to
    // carry out and realign
    if(i32<CL) {
      uint32_t a = CT[i32];
      uint32_t c = a - carry;
      carry = (c>a);
      CT[i32] = c;
    }
    i64++; // We have calculated another 64 bits.
  }

  // Now carry out
  for(;i64<CL64;i64++) {
    uint64_t a = CT64[i64];
    uint64_t c = a - carry;
    carry = (c>a);
    CT64[i64] = c;

    // If there is nothing more to carry, we can return early;
    if(carry == 0) {
      return;
    }
  }

  // Now another 32 bits if there is still some carry left
  if(CL%2 == 1) {
    i32 = i64*2;

    uint32_t a = CT[i32];
    uint32_t c = a - carry;
    carry = (c>a);
    CT[i32] = c;
    i32++;
  }

  // if we are at the end and there is still carry, we fucked up.
  if(carry) {
    printf("Carry Error\n");
    abort();
  }

  // FIXME: handle CL%2 == 1

  // Original code
  /*
  int32_t carry = 0;
  for (size_t i = 0; i<CL; i++) {
    int64_t word = (int64_t)CT[i]-(int64_t)carry;
    if(i<AL) {
      word -= (int64_t)AT[i];
    }
    carry = 0;
    if (word < 0){
      word += 0x100000000;
      carry = 1;
    }
    CT[i] = word;
  }
  */
}

#define SIZE1 (10001)
#define SIZE2 (8001)

int main(void) {
  size_t CL = SIZE1;
  uint32_t* CT = malloc(CL*sizeof(uint32_t));
  for(int i=0;i<CL;i++) {
    CT[i] = 0xffff+i;
  }

  size_t AL = SIZE2;
  uint32_t* AT = malloc(AL*sizeof(uint32_t));
  for(int i=0;i<AL;i++) {
    AT[i] = 500-i;
  }
  AT[0] = 0;

  for(int i=0;i<100;i++) {
    printf("Iter %d\n", i);
    uint32_t ccb = checksum(CT, CL);
    uint32_t ca  = checksum(AT, AL);
    _ll_sub_inplace(CT, CL, AT, AL);
    if(checksum_sub(ccb, ca) != checksum(CT, CL)) {
      printf("Checksum error\n");
      abort();
    }
  }

  printf("%08x\n", CT[7]); // prevent dead code removal
}

