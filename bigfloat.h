#define YCL_BIGFLOAT_EXTRA_PRECISION    2

typedef struct __BigFloat {
    bool sign;      //  true = positive or zero, false = negative
    int64_t exp;    //  Exponent
    size_t L;       //  Length
    uint32_t* T;
} BigFloat;

void bigfloat_negate(BigFloat& num);
void bigfloat_set(BigFloat &target, uint32_t x, bool sign_);
uint32_t _bigfloat_word_at(const BigFloat &target, int64_t mag);
int _bigfloat_ucmp(const BigFloat &a, const BigFloat &b);
void bigfloat_add(BigFloat &target, const BigFloat &a, const BigFloat &b, size_t p);
void bigfloat_sub(BigFloat &target, const BigFloat &a, const BigFloat &b, size_t p);
void bigfloat_mul(BigFloat &target, const BigFloat &a, const BigFloat &b, size_t p, int tds);
void bigfloat_div(BigFloat &target, const BigFloat &a, const BigFloat &b, size_t p, int tds);
void bigfloat_rcp(BigFloat &target, const BigFloat &a, size_t p, int tds);
std::string bigfloat_to_string(const BigFloat& value, size_t digits);
void bigfloat_free(BigFloat &target);
void bigfloat_new(BigFloat &target);
