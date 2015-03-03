#define YCL_BIGFLOAT_EXTRA_PRECISION    2
class BigFloat{
public:
    BigFloat(BigFloat &&x);
    BigFloat& operator=(BigFloat &&x);
    ~BigFloat();

    BigFloat();

    std::string to_string    (size_t digits = 0) const;
    std::string to_string_sci(size_t digits = 0) const;
    uint32_t word_at(int64_t mag) const;

    BigFloat add(const BigFloat &x,size_t p = 0) const;
    BigFloat sub(const BigFloat &x,size_t p = 0) const;
    BigFloat mul(const BigFloat &x,size_t p = 0,int threads = 0) const;
    BigFloat rcp(size_t p,int threads = 0) const;
    BigFloat div(const BigFloat &x,size_t p,int threads = 0) const;

    bool sign;      //  true = positive or zero, false = negative
    int64_t exp;    //  Exponent
    size_t L;       //  Length
    uint32_t* T;

private:
    //  Internal helpers
    int64_t to_string_trimmed(size_t digits,std::string &str) const;
    BigFloat uadd(const BigFloat &x,size_t p) const;
    BigFloat usub(const BigFloat &x,size_t p) const;

    friend BigFloat invsqrt(uint32_t x,size_t p,int threads);
};
BigFloat invsqrt(uint32_t x,size_t p);

void bigfloat_negate(BigFloat& num);
void bigfloat_set(BigFloat &target, uint32_t x, bool sign_);
uint32_t _bigfloat_word_at(const BigFloat &target, int64_t mag);
int _bigfloat_ucmp(const BigFloat &a, const BigFloat &b);
