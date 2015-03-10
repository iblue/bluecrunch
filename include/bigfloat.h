#include <stdint.h>

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

// Arithmatic
void bigfloat_neg(bigfloat_t target);
void bigfloat_add(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_sub(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_div(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_rcp(bigfloat_t target, const bigfloat_t a,                     size_t p, int tds);

// Output
size_t bigfloat_to_string(char* string, const bigfloat_t value, size_t digits);
void bigfloat_print(const char* name,   const bigfloat_t value);
