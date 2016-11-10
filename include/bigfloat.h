#include <stdint.h>
#include <stddef.h>

#define YCL_BIGFLOAT_EXTRA_PRECISION    2

typedef struct {
    uint32_t* coef;    // Coefficients
    int64_t   exp;     // Exponent
    size_t    len;     // Length
    int       sign;    // 1: pos, 0: neg
} _bigfloat_struct;

typedef _bigfloat_struct bigfloat_t[1];


// Bigfloat management and setting
void bigfloat_zero  (      bigfloat_t target);
int  bigfloat_iszero(const bigfloat_t target);
void bigfloat_set   (      bigfloat_t target, uint32_t value);
void bigfloat_free  (      bigfloat_t target);
void bigfloat_new   (      bigfloat_t target);
void bigfloat_copy  (      bigfloat_t target, const bigfloat_t source);
void bigfloat_alloc(bigfloat_t target, size_t size);
void bigfloat_realloc(bigfloat_t target, size_t size);

// Arithmatic
void bigfloat_neg(bigfloat_t target);
void bigfloat_floor(bigfloat_t target);
void bigfloat_add(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_sub(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_mului(bigfloat_t a, uint32_t b);
void bigfloat_div(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_rcp(bigfloat_t target, const bigfloat_t a,                     size_t p, int tds);
void bigfloat_exp(bigfloat_t target, const bigfloat_t a, const uint64_t exp, int tds);
size_t bigfloat_radix_decimals(bigfloat_t target);

// Radix conversion
void bigfloat_radix(bigfloat_t target, const bigfloat_t a, int tds);

// Output
size_t bigfloat_to_string(char* string, const bigfloat_t value, size_t digits, int base);
void bigfloat_print(const char* name,   const bigfloat_t value);
void bigfloat_print10(const char* name, const bigfloat_t value);
