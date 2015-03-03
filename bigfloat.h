#define YCL_BIGFLOAT_EXTRA_PRECISION    2
class BigFloat{
public:
    BigFloat(BigFloat &&x);
    BigFloat& operator=(BigFloat &&x);
    ~BigFloat();

    BigFloat();
    BigFloat(uint32_t x,bool sign = true);

    std::string to_string    (size_t digits = 0) const;
    std::string to_string_sci(size_t digits = 0) const;
    uint32_t word_at(int64_t mag) const;

    BigFloat mul(uint32_t x) const;
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
    int ucmp(const BigFloat &x) const;
    BigFloat uadd(const BigFloat &x,size_t p) const;
    BigFloat usub(const BigFloat &x,size_t p) const;

    friend BigFloat invsqrt(uint32_t x,size_t p,int threads);
};
BigFloat invsqrt(uint32_t x,size_t p);

void bigfloat_negate(BigFloat& num);
