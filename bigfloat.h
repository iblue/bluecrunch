#define YCL_BIGFLOAT_EXTRA_PRECISION    2

typedef struct {
    int sign;       //  1 = positive or zero, 0 = negative
    int64_t exp;    //  Exponent
    size_t L;       //  Length
    uint32_t* T;
} _bigfloat_struct;

typedef _bigfloat_struct bigfloat_t[1];

void bigfloat_negate(bigfloat_t num);
void bigfloat_set(bigfloat_t target, uint32_t x, int sign_);
uint32_t _bigfloat_word_at(const bigfloat_t target, int64_t mag);
int _bigfloat_ucmp(const bigfloat_t a, const bigfloat_t b);
void bigfloat_add(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_sub(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p);
void bigfloat_mul(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_div(bigfloat_t target, const bigfloat_t a, const bigfloat_t b, size_t p, int tds);
void bigfloat_rcp(bigfloat_t target, const bigfloat_t a, size_t p, int tds);
size_t bigfloat_to_string(char* string, const bigfloat_t value, size_t digits);
void bigfloat_free(bigfloat_t target);
void bigfloat_new(bigfloat_t target);
