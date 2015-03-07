#define YCL_BIGFLOAT_EXTRA_PRECISION    2

typedef struct {
    int sign;       //  1 = positive or zero, 0 = negative
    int64_t exp;    //  Exponent
    size_t L;       //  Length
    uint32_t* T;
} _bigfloat_struct;

typedef _bigfloat_struct BigFloat[1];

void bigfloat_negate(BigFloat num);
void bigfloat_set(BigFloat target, uint32_t x, int sign_);
uint32_t _bigfloat_word_at(const BigFloat target, int64_t mag);
int _bigfloat_ucmp(const BigFloat a, const BigFloat b);
void bigfloat_add(BigFloat target, const BigFloat a, const BigFloat b, size_t p);
void bigfloat_sub(BigFloat target, const BigFloat a, const BigFloat b, size_t p);
void bigfloat_mul(BigFloat target, const BigFloat a, const BigFloat b, size_t p, int tds);
void bigfloat_div(BigFloat target, const BigFloat a, const BigFloat b, size_t p, int tds);
void bigfloat_rcp(BigFloat target, const BigFloat a, size_t p, int tds);
size_t bigfloat_to_string(char* string, const BigFloat value, size_t digits);
void bigfloat_free(BigFloat target);
void bigfloat_new(BigFloat target);
