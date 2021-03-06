#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <pmmintrin.h>
#include <immintrin.h> // More Magic!
#include <cilk/cilk.h>
#include "fft.h"
#include "intrinsic.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define LENGTHS 125

size_t lengths[LENGTHS] = {
  2,                   // 0 => 2*2^0
  3,                   // 1 => 3*2^0
  4,                   // 2 => 2*2^1
  6,                   // 3 => 3*2^1
  8,                   // 4 => 2*2^2
  12,                  // 5 => 3*2^2
  16,                  // 6 => 2*2^3
  24,                  // 7 => 3*2^3
  32,                  // 8 => 2*2^4
  48,                  // 9 => 3*2^4
  64,                  // 10 => 2*2^5
  96,                  // 11 => 3*2^5
  128,                 // 12 => 2*2^6
  192,                 // 13 => 3*2^6
  256,                 // 14 => 2*2^7
  384,                 // 15 => 3*2^7
  512,                 // 16 => 2*2^8
  768,                 // 17 => 3*2^8
  1024,                // 18 => 2*2^9
  1536,                // 19 => 3*2^9
  2048,                // 20 => 2*2^10
  3072,                // 21 => 3*2^10
  4096,                // 22 => 2*2^11
  6144,                // 23 => 3*2^11
  8192,                // 24 => 2*2^12
  12288,               // 25 => 3*2^12
  16384,               // 26 => 2*2^13
  24576,               // 27 => 3*2^13
  32768,               // 28 => 2*2^14
  49152,               // 29 => 3*2^14
  65536,               // 30 => 2*2^15
  98304,               // 31 => 3*2^15
  131072,              // 32 => 2*2^16
  196608,              // 33 => 3*2^16
  262144,              // 34 => 2*2^17
  393216,              // 35 => 3*2^17
  524288,              // 36 => 2*2^18
  786432,              // 37 => 3*2^18
  1048576,             // 38 => 2*2^19
  1572864,             // 39 => 3*2^19
  2097152,             // 40 => 2*2^20
  3145728,             // 41 => 3*2^20
  4194304,             // 42 => 2*2^21
  6291456,             // 43 => 3*2^21
  8388608,             // 44 => 2*2^22
  12582912,            // 45 => 3*2^22
  16777216,            // 46 => 2*2^23
  25165824,            // 47 => 3*2^23
  33554432,            // 48 => 2*2^24
  50331648,            // 49 => 3*2^24
  67108864,            // 50 => 2*2^25
  100663296,           // 51 => 3*2^25
  134217728,           // 52 => 2*2^26
  201326592,           // 53 => 3*2^26
  268435456,           // 54 => 2*2^27
  402653184,           // 55 => 3*2^27
  536870912,           // 56 => 2*2^28
  805306368,           // 57 => 3*2^28
  1073741824,          // 58 => 2*2^29
  1610612736,          // 59 => 3*2^29
  2147483648,          // 60 => 2*2^30
  3221225472,          // 61 => 3*2^30
  4294967296,          // 62 => 2*2^31
  6442450944,          // 63 => 3*2^31
  8589934592,          // 64 => 2*2^32
  12884901888,         // 65 => 3*2^32
  17179869184,         // 66 => 2*2^33
  25769803776,         // 67 => 3*2^33
  34359738368,         // 68 => 2*2^34
  51539607552,         // 69 => 3*2^34
  68719476736,         // 70 => 2*2^35
  103079215104,        // 71 => 3*2^35
  137438953472,        // 72 => 2*2^36
  206158430208,        // 73 => 3*2^36
  274877906944,        // 74 => 2*2^37
  412316860416,        // 75 => 3*2^37
  549755813888,        // 76 => 2*2^38
  824633720832,        // 77 => 3*2^38
  1099511627776,       // 78 => 2*2^39
  1649267441664,       // 79 => 3*2^39
  2199023255552,       // 80 => 2*2^40
  3298534883328,       // 81 => 3*2^40
  4398046511104,       // 82 => 2*2^41
  6597069766656,       // 83 => 3*2^41
  8796093022208,       // 84 => 2*2^42
  13194139533312,      // 85 => 3*2^42
  17592186044416,      // 86 => 2*2^43
  26388279066624,      // 87 => 3*2^43
  35184372088832,      // 88 => 2*2^44
  52776558133248,      // 89 => 3*2^44
  70368744177664,      // 90 => 2*2^45
  105553116266496,     // 91 => 3*2^45
  140737488355328,     // 92 => 2*2^46
  211106232532992,     // 93 => 3*2^46
  281474976710656,     // 94 => 2*2^47
  422212465065984,     // 95 => 3*2^47
  562949953421312,     // 96 => 2*2^48
  844424930131968,     // 97 => 3*2^48
  1125899906842624,    // 98 => 2*2^49
  1688849860263936,    // 99 => 3*2^49
  2251799813685248,    // 100 => 2*2^50
  3377699720527872,    // 101 => 3*2^50
  4503599627370496,    // 102 => 2*2^51
  6755399441055744,    // 103 => 3*2^51
  6755399441055744,    // 103 => 3*2^51
  9007199254740992,    // 104 => 2*2^52
  13510798882111488,   // 105 => 3*2^52
  18014398509481984,   // 106 => 2*2^53
  27021597764222976,   // 107 => 3*2^53
  36028797018963968,   // 108 => 2*2^54
  54043195528445952,   // 109 => 3*2^54
  72057594037927936,   // 110 => 2*2^55
  108086391056891904,  // 111 => 3*2^55
  144115188075855872,  // 112 => 2*2^56
  216172782113783808,  // 113 => 3*2^56
  288230376151711744,  // 114 => 2*2^57
  432345564227567616,  // 115 => 3*2^57
  576460752303423488,  // 116 => 2*2^58
  864691128455135232,  // 117 => 3*2^58
  1152921504606846976, // 118 => 2*2^59
  1729382256910270464, // 119 => 3*2^59
  2305843009213693952, // 120 => 2*2^60
  3458764513820540928, // 121 => 3*2^60
  4611686018427387904, // 122 => 2*2^61
  6917529027641081856, // 123 => 3*2^61
  /*
  9223372036854775808, // 124 => 2*2^62
  13835058055282163712 // 125 => 3*2^62
  */
};

complex double* twiddle_table[LENGTHS];
int twiddle_table_size = 0;

size_t _table_select(size_t length) {
  size_t p     = __builtin_popcountl(length); // 1 for 2*2^k, 2 for 3*2^k
  size_t k     = __builtin_ctzl(length); // returns k or k+1 for 2*2^k
#ifdef DEBUG
  p--; // 0 for 2*2^k, 1 for 3*2^k
  k -= p ? 0 : 1; // correct k by -1 if 2*2^k
  size_t base  = 1 << k; // 2^k
  size_t ret   = base << 1; // 2*2^k
  ret += p ? base : 0; // adds 2^k if p == 1 (leading 2*2^k if p == 0, 3*2^k else)

  // did we get it right? => return
  if(ret == length) {
    return (k<<1)+p;
  } else {
    fprintf(stderr, "Invalid table select for len %ld\n", length);
    abort();
  }
#else
  if(p==1) {
    return 2*(k-1);
  } else {
    return 2*k+1;
  }
#endif
}

size_t table_select(size_t length) {
  size_t index = _table_select(length);
  if(index > twiddle_table_size) {
    return -1;
  }
  return index;
}

void _build_table(size_t length) {
  double omega = 2 * M_PI / length;
  length /= 2;

  //  Build the sub-table.
  complex double *sub_table = (complex double*)_mm_malloc(length*sizeof(complex double), 32);

  cilk_for(size_t c = 0; c < length; c++) {
    //  Generate Twiddle Factor
    double angle = omega * c;
    complex double twiddle_factor = cos(angle)+I*sin(angle);
    sub_table[c] = twiddle_factor;
  }

#ifdef DEBUG
  printf("Pushing table length %ld to index %ld\n", length, _table_select(length));
#endif

  //  Push into main table.
  twiddle_table[_table_select(length)] = sub_table;
}

void fft_ensure_table(size_t length) {
  size_t k=2;
  while(1) {
    if(2*k > length) {
      break;
    }
    cilk_spawn _build_table(2*k);
    if(3*k > length) {
      break;
    }
    cilk_spawn _build_table(3*k);
    k*=2;
  }

  cilk_sync;

  twiddle_table_size = _table_select(k)+1;
}

void fft_free_table() {
  for(size_t i=0;i<twiddle_table_size;i++) {
    _mm_free(twiddle_table[i]);
  }
}
