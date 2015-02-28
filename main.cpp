#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
using std::cout;
using std::endl;

//  SIMD
#include <malloc.h>
#include <pmmintrin.h>

#include <omp.h>

#if USE_CHRONO
#include <chrono>
#else
#include <time.h>
#endif

#include "fft.h"

namespace Bluecrunch {
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Helpers
double wall_clock(){
    //  Get the clock in seconds.
#if USE_CHRONO
    auto ratio_object = std::chrono::high_resolution_clock::period();
    double ratio = (double)ratio_object.num / ratio_object.den;
    return std::chrono::high_resolution_clock::now().time_since_epoch().count() * ratio;
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}
void dump_to_file(const char *path,const std::string &str){
    //  Dump a string to a file.

    FILE *file = fopen(path,"wb");
    if (file == NULL)
        throw "Cannot Create File";

    fwrite(str.c_str(),1,str.size(),file);
    fclose(file);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  BigFloat object
/*  This is the big floating-point object. It represents an arbitrary precision
 *  floating-point number.
 * 
 *  Its numerical value is equal to:
 * 
 *      word = 10^9
 *      word^exp * (T[0] + T[1]*word + T[2]*word^2 + ... + T[L - 1]*word^(L - 1))
 * 
 *  T is an array of 32-bit integers. Each integer stores 9 decimal digits
 *  and must always have a value in the range [0, 999999999].
 * 
 *  T[L - 1] must never be zero. 
 * 
 *  The number is positive when (sign = true) and negative when (sign = false).
 *  Zero is represented as (sign = true) and (L = 0).
 * 
 */
#define YCL_BIGFLOAT_EXTRA_PRECISION    2
class BigFloat{
public:
    BigFloat(BigFloat &&x);
    BigFloat& operator=(BigFloat &&x);

    BigFloat();
    BigFloat(uint32_t x,bool sign = true);

    std::string to_string    (size_t digits = 0) const;
    std::string to_string_sci(size_t digits = 0) const;
    size_t get_precision() const;
    int64_t get_exponent() const;
    uint32_t word_at(int64_t mag) const;

    void negate();
    BigFloat mul(uint32_t x) const;
    BigFloat add(const BigFloat &x,size_t p = 0) const;
    BigFloat sub(const BigFloat &x,size_t p = 0) const;
    BigFloat mul(const BigFloat &x,size_t p = 0,int threads = 0) const;
    BigFloat rcp(size_t p,int threads = 0) const;
    BigFloat div(const BigFloat &x,size_t p,int threads = 0) const;

private:
    bool sign;      //  true = positive or zero, false = negative
    int64_t exp;    //  Exponent
    size_t L;       //  Length
    std::unique_ptr<uint32_t[]> T;

    //  Internal helpers
    int64_t to_string_trimmed(size_t digits,std::string &str) const;
    int ucmp(const BigFloat &x) const;
    BigFloat uadd(const BigFloat &x,size_t p) const;
    BigFloat usub(const BigFloat &x,size_t p) const;

    friend BigFloat invsqrt(uint32_t x,size_t p,int threads);
};
BigFloat invsqrt(uint32_t x,size_t p);
////////////////////////////////////////////////////////////////////////////////
//  Move operators
BigFloat::BigFloat(BigFloat &&x)
    : sign(x.sign)
    , exp(x.exp)
    , L(x.L)
    , T(std::move(x.T))
{
    x.sign  = true;
    x.exp   = 0;
    x.L     = 0;
}
BigFloat& BigFloat::operator=(BigFloat &&x){
    sign    = x.sign;
    exp     = x.exp;
    L       = x.L;
    T       = std::move(x.T);

    x.sign  = true;
    x.exp   = 0;
    x.L     = 0;

    return *this;
}
////////////////////////////////////////////////////////////////////////////////
//  Constructors
BigFloat::BigFloat()
    : sign(true)
    , exp(0)
    , L(0)
{}
BigFloat::BigFloat(uint32_t x,bool sign_)
    : sign(true)
    , exp(0)
    , L(1)
{
    //  Construct a BigFloat with a value of x and the specified sign.

    if (x == 0){
        L = 0;
        return;
    }
    sign = sign_;

    T = std::unique_ptr<uint32_t[]>(new uint32_t[1]);
    T[0] = x;
}
////////////////////////////////////////////////////////////////////////////////
//  String Conversion
int64_t BigFloat::to_string_trimmed(size_t digits,std::string &str) const{
    //  Converts this object to a string with "digits" significant figures.

    //  After calling this function, the following expression is equal to the
    //  numeric value of this object. (after truncation of precision)
    //      str + " * 10^" + (return value)

    if (L == 0){
        str = "0";
        return 0;
    }

    //  Collect operands
    int64_t exponent = exp;
    size_t length = L;
    uint32_t *ptr = T.get();

    if (digits == 0){
        //  Use all digits.
        digits = length * 9;
    }else{
        //  Truncate precision
        size_t words = (digits + 17) / 9;
        if (words < length){
            size_t chop = length - words;
            exponent += chop;
            length = words;
            ptr += chop;
        }
    }
    exponent *= 9;

    //  Build string
    char buffer[] = "012345678";
    str.clear();
    size_t c = length;
    while (c-- > 0){
        uint32_t word = ptr[c];
        for (int i = 8; i >= 0; i--){
            buffer[i] = word % 10 + '0';
            word /= 10;
        }
        str += buffer;
    }

    //  Count leading zeros
    size_t leading_zeros = 0;
    while (str[leading_zeros] == '0')
        leading_zeros++;
    digits += leading_zeros;

    //  Truncate
    if (digits < str.size()){
        exponent += str.size() - digits;
        str.resize(digits);
    }

    return exponent;
}
std::string BigFloat::to_string(size_t digits) const{
    //  Convert this number to a string. Auto-select format type.
    if (L == 0)
        return "0.";

    int64_t mag = exp + L;

    //  Use scientific notation of out of range.
    if (mag > 1 || mag < 0)
        return to_string_sci();

    //  Convert
    std::string str;
    int64_t exponent = to_string_trimmed(digits,str);

    //  Less than 1
    if (mag == 0){
        if (sign)
            return std::string("0.") + str;
        else
            return std::string("-0.") + str;
    }

    //  Get a string with the digits before the decimal place.
    std::string before_decimal = std::to_string((long long)T[L - 1]);

    //  Nothing after the decimal place.
    if (exponent >= 0){
        if (sign){
            return before_decimal + ".";
        }else{
            return std::string("-") + before_decimal + ".";
        }
    }

    //  Get digits after the decimal place.
    std::string after_decimal = str.substr((size_t)(str.size() + exponent),(size_t)-exponent);

    if (sign){
        return before_decimal + "." + after_decimal;
    }else{
        return std::string("-") + before_decimal + "." + after_decimal;
    }
}
std::string BigFloat::to_string_sci(size_t digits) const{
    //  Convert to string in scientific notation.
    if (L == 0)
        return "0.";

    //  Convert
    std::string str;
    int64_t exponent = to_string_trimmed(digits,str);

    //  Strip leading zeros.
    {
        size_t leading_zeros = 0;
        while (str[leading_zeros] == '0')
            leading_zeros++;
        str = &str[leading_zeros];
    }

    //  Insert decimal place
    exponent += str.size() - 1;
    str = str.substr(0,1) + "." + &str[1];

    //  Add exponent
    if (exponent != 0){
        str += " * 10^";
        str += std::to_string(exponent);
    }

    //  Add sign
    if (!sign)
        str = std::string("-") + str;

    return str;
}
////////////////////////////////////////////////////////////////////////////////
//  Getters
size_t BigFloat::get_precision() const{
    //  Returns the precision of the number in words.
    //  Note that each word is 9 decimal digits.
    return L;
}
int64_t BigFloat::get_exponent() const{
    //  Returns the exponent of the number in words.
    //  Note that each word is 9 decimal digits.
    return exp;
}
uint32_t BigFloat::word_at(int64_t mag) const{
    //  Returns the word at the mag'th digit place.
    //  This is useful for additions where you need to access a specific "digit place"
    //  of the operand without having to worry if it's out-of-bounds.

    //  This function is mathematically equal to:
    //      (return value) = floor(this * (10^9)^-mag) % 10^9

    if (mag < exp)
        return 0;
    if (mag >= exp + (int64_t)L)
        return 0;
    return T[(size_t)(mag - exp)];
}
int BigFloat::ucmp(const BigFloat &x) const{
    //  Compare function that ignores the sign.
    //  This is needed to determine which direction subtractions will go.

    //  Magnitude
    int64_t magA = exp + L;
    int64_t magB = x.exp + x.L;
    if (magA > magB)
        return 1;
    if (magA < magB)
        return -1;

    //  Compare
    int64_t mag = magA;
    while (mag >= exp || mag >= x.exp){
        uint32_t wordA = word_at(mag);
        uint32_t wordB = x.word_at(mag);
        if (wordA < wordB)
            return -1;
        if (wordA > wordB)
            return 1;
        mag--;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
//  Arithmetic
void BigFloat::negate(){
    //  Negate this number.
    if (L == 0)
        return;
    sign = !sign;
}
BigFloat BigFloat::mul(uint32_t x) const{
    //  Multiply by a 32-bit unsigned integer.

    if (L == 0 || x == 0)
        return BigFloat();

    //  Compute basic fields.
    BigFloat z;
    z.sign = sign;
    z.exp  = exp;
    z.L    = L;

    //  Allocate mantissa
    z.T = std::unique_ptr<uint32_t[]>(new uint32_t[z.L + 1]);

    uint64_t carry = 0;
    for (size_t c = 0; c < L; c++){
        carry += (uint64_t)T[c] * x;                //  Multiply and add to carry
        z.T[c] = (uint32_t)(carry % 1000000000);    //  Store bottom 9 digits
        carry /= 1000000000;                        //  Shift down the carry
    }

    //  Carry out
    if (carry != 0)
        z.T[z.L++] = (uint32_t)carry;

    return z;
}
BigFloat BigFloat::uadd(const BigFloat &x,size_t p) const{
    //  Perform addition ignoring the sign of the two operands.

    //  Magnitude
    int64_t magA = exp + L;
    int64_t magB = x.exp + x.L;
    int64_t top = std::max(magA,magB);
    int64_t bot = std::min(exp,x.exp);

    //  Target length
    int64_t TL = top - bot;

    if (p == 0){
        //  Default value. No trunction.
        p = (size_t)TL;
    }else{
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    //  Perform precision truncation.
    if (TL > (int64_t)p){
        bot = top - p;
        TL = p;
    }

    //  Compute basic fields.
    BigFloat z;
    z.sign  = sign;
    z.exp   = bot;
    z.L     = (uint32_t)TL;

    //  Allocate mantissa
    z.T = std::unique_ptr<uint32_t[]>(new uint32_t[z.L + 1]);

    //  Add
    uint32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++){
        uint32_t word = word_at(bot) + x.word_at(bot) + carry;
        carry = 0;
        if (word >= 1000000000){
            word -= 1000000000;
            carry = 1;
        }
        z.T[c] = word;
    }

    //  Carry out
    if (carry != 0){
        z.T[z.L++] = 1;
    }

    return z;
}
BigFloat BigFloat::usub(const BigFloat &x,size_t p) const{
    //  Perform subtraction ignoring the sign of the two operands.

    //  "this" must be greater than or equal to x. Otherwise, the behavior
    //  is undefined.

    //  Magnitude
    int64_t magA = exp + L;
    int64_t magB = x.exp + x.L;
    int64_t top = std::max(magA,magB);
    int64_t bot = std::min(exp,x.exp);

    //  Truncate precision
    int64_t TL = top - bot;

    if (p == 0){
        //  Default value. No trunction.
        p = (size_t)TL;
    }else{
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    if (TL > (int64_t)p){
        bot = top - p;
        TL = p;
    }

    //  Compute basic fields.
    BigFloat z;
    z.sign  = sign;
    z.exp   = bot;
    z.L     = (uint32_t)TL;

    //  Allocate mantissa
    z.T = std::unique_ptr<uint32_t[]>(new uint32_t[z.L]);

    //  Subtract
    int32_t carry = 0;
    for (size_t c = 0; bot < top; bot++, c++){
        int32_t word = (int32_t)word_at(bot) - (int32_t)x.word_at(bot) - carry;
        carry = 0;
        if (word < 0){
            word += 1000000000;
            carry = 1;
        }
        z.T[c] = word;
    }

    //  Strip leading zeros
    while (z.L > 0 && z.T[z.L - 1] == 0)
        z.L--;
    if (z.L == 0){
        z.exp = 0;
        z.sign = true;
        z.T.reset();
    }

    return z;
}
BigFloat BigFloat::add(const BigFloat &x,size_t p) const{
    //  Addition

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Same sign. Add.
    if (sign == x.sign)
        return uadd(x,p);

    //  this > x
    if (ucmp(x) > 0)
        return usub(x,p);

    //  this < x
    BigFloat z = x.usub(*this,p);
    z.negate();
    return z;
}
BigFloat BigFloat::sub(const BigFloat &x,size_t p) const{
    //  Subtraction

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Different sign. Add.
    if (sign != x.sign)
        return uadd(x,p);

    //  this > x
    if (ucmp(x) > 0)
        return usub(x,p);

    //  this < x
    BigFloat z = x.usub(*this,p);
    z.negate();
    return z;
}
BigFloat BigFloat::mul(const BigFloat &x,size_t p,int tds) const{
    //  Multiplication

    //  The target precision is p.
    //  If (p = 0), then no truncation is done. The entire operation is done
    //  at maximum precision with no data loss.

    //  Either operand is zero.
    if (L == 0 || x.L == 0)
        return BigFloat();

    if (p == 0){
        //  Default value. No trunction.
        p = L + x.L;
    }else{
        //  Increase precision
        p += YCL_BIGFLOAT_EXTRA_PRECISION;
    }

    //  Collect operands.
    int64_t Aexp = exp;
    int64_t Bexp = x.exp;
    size_t AL = L;
    size_t BL = x.L;
    uint32_t *AT = T.get();
    uint32_t *BT = x.T.get();

    //  Perform precision truncation.
    if (AL > p){
        size_t chop = AL - p;
        AL = p;
        Aexp += chop;
        AT += chop;
    }
    if (BL > p){
        size_t chop = BL - p;
        BL = p;
        Bexp += chop;
        BT += chop;
    }

    //  Compute basic fields.
    BigFloat z;
    z.sign = sign == z.sign;    //  Sign is positive if signs are equal.
    z.exp  = Aexp + Bexp;       //  Add the exponents.
    z.L    = AL + BL;           //  Add the lenghts for now. May need to correct later.

    //  Allocate mantissa
    z.T = std::unique_ptr<uint32_t[]>(new uint32_t[z.L]);

    //  Perform multiplication.

    //  Determine minimum FFT size.
    int k = 0;
    size_t length = 1;
    while (length < 3*z.L){
        length <<= 1;
        k++;
    }

    //  Perform a convolution using FFT.
    //  Yeah, this is slow for small sizes, but it's asympotically optimal.

    //  3 digits per point is small enough to not encounter round-off error
    //  until a transform size of 2^30.
    //  A transform length of 2^29 allows for the maximum product size to be
    //  2^29 * 3 = 1,610,612,736 decimal digits.
    if (k > 29)
        throw "FFT size limit exceeded.";

    //  Allocate FFT arrays
    SIMD_delete deletor;
    auto Ta = std::unique_ptr<__m128d[],SIMD_delete>((__m128d*)_mm_malloc(length * sizeof(__m128d),16),deletor);
    auto Tb = std::unique_ptr<__m128d[],SIMD_delete>((__m128d*)_mm_malloc(length * sizeof(__m128d),16),deletor);

    //  Make sure the twiddle table is big enough.
    if ((int)twiddle_table.size() - 1 < k)
        throw "Table is not large enough.";

    int_to_fft(Ta.get(),k,AT,AL);           //  Convert 1st operand
    int_to_fft(Tb.get(),k,BT,BL);           //  Convert 2nd operand
    fft_forward(Ta.get(),k,tds);            //  Transform 1st operand
    fft_forward(Tb.get(),k,tds);            //  Transform 2nd operand
    fft_pointwise(Ta.get(),Tb.get(),k);     //  Pointwise multiply
    fft_inverse(Ta.get(),k);                //  Perform inverse transform.
    fft_to_int(Ta.get(),k,z.T.get(),z.L);   //  Convert back to word array.

    //  Check top word and correct length.
    if (z.T[z.L - 1] == 0)
        z.L--;

    return z;
}
BigFloat BigFloat::rcp(size_t p,int tds) const{
    //  Compute reciprocal using Newton's Method.

    //  r1 = r0 - (r0 * x - 1) * r0

    if (L == 0)
        throw "Divide by Zero";

    //  Collect operand
    int64_t Aexp = exp;
    size_t AL = L;
    uint32_t *AT = T.get();

    //  End of recursion. Generate starting point.
    if (p == 0){
        //  Truncate precision to 3.
        p = 3;
        if (AL > p){
            size_t chop = AL - p;
            AL = p;
            Aexp += chop;
            AT += chop;
        }

        //  Convert number to floating-point.
        double val = AT[0];
        if (AL >= 2)
            val += AT[1] * 1000000000.;
        if (AL >= 3)
            val += AT[2] * 1000000000000000000.;

        //  Compute reciprocal.
        val = 1. / val;
        Aexp = -Aexp;

        //  Scale
        while (val < 1000000000.){
            val *= 1000000000.;
            Aexp--;
        }

        //  Rebuild a BigFloat.
        uint64_t val64 = (uint64_t)val;

        BigFloat out;
        out.sign = sign;

        out.T = std::unique_ptr<uint32_t[]>(new uint32_t[2]);
        out.T[0] = (uint32_t)(val64 % 1000000000);
        out.T[1] = (uint32_t)(val64 / 1000000000);
        out.L = 2;
        out.exp = Aexp;

        return out;
    }

    //  Half the precision
    size_t s = p / 2 + 1;
    if (p == 1) s = 0;
    if (p == 2) s = 1;

    //  Recurse at half the precision
    BigFloat T = rcp(s,tds);

    //  r1 = r0 - (r0 * x - 1) * r0
    return T.sub(this->mul(T,p,tds).sub(BigFloat(1),p).mul(T,p,tds),p);
}
BigFloat BigFloat::div(const BigFloat &x,size_t p,int tds) const{
    //  Division
    return this->mul(x.rcp(p,tds),p,tds);
}
BigFloat invsqrt(uint32_t x,size_t p,int tds){
    //  Compute inverse square root using Newton's Method.

    //            (  r0^2 * x - 1  )
    //  r1 = r0 - (----------------) * r0
    //            (       2        )

    if (x == 0)
        throw "Divide by Zero";

    //  End of recursion. Generate starting point.
    if (p == 0){
        double val = 1. / sqrt((double)x);

        int64_t exponent = 0;

        //  Scale
        while (val < 1000000000.){
            val *= 1000000000.;
            exponent--;
        }

        //  Rebuild a BigFloat.
        uint64_t val64 = (uint64_t)val;

        BigFloat out;
        out.sign = true;

        out.T = std::unique_ptr<uint32_t[]>(new uint32_t[2]);
        out.T[0] = (uint32_t)(val64 % 1000000000);
        out.T[1] = (uint32_t)(val64 / 1000000000);
        out.L = 2;
        out.exp = exponent;

        return out;
    }

    //  Half the precision
    size_t s = p / 2 + 1;
    if (p == 1) s = 0;
    if (p == 2) s = 1;

    //  Recurse at half the precision
    BigFloat T = invsqrt(x,s,tds);

    BigFloat temp = T.mul(T,p);     //  r0^2
    temp = temp.mul(x,p,tds);       //  r0^2 * x
    temp = temp.sub(BigFloat(1),p); //  r0^2 * x - 1
    temp = temp.mul(500000000);         //  (r0^2 * x - 1) / 2
    temp.exp--;
    temp = temp.mul(T,p,tds);       //  (r0^2 * x - 1) / 2 * r0
    return T.sub(temp,p);           //  r0 - (r0^2 * x - 1) / 2 * r0
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  e
double logf_approx(double x){
    //  Returns a very good approximation to log(x!).
    //  log(x!) ~ (x + 1/2) * (log(x) - 1) + (log(2*pi) + 1) / 2
    //  This approximation gets better as x is larger.
    if (x <= 1) return 0;
    return (x + .5) * (log(x) - 1) + (1.4189385332046727417803297364056176398613974736378);
}
size_t e_terms(size_t p){
    //  Returns the # of terms needed to reach a precision of p.

    //  The taylor series converges to log(x!) / log(10) decimal digits after
    //  x terms. So to find the number of terms needed to reach a precision of p
    //  we need to solve this question for x:
    //      p = log(x!) / log(1000000000)

    //  This function solves this equation via binary search.

    double sizeL = (double)p * 20.723265836946411156161923092159277868409913397659 + 1;

    size_t a = 0;
    size_t b = 1;

    //  Double up
    while (logf_approx((double)b) < sizeL)
        b <<= 1;

    //  Binary search
    while (b - a > 1){
        size_t m = (a + b) >> 1;

        if (logf_approx((double)m) < sizeL)
            a = m;
        else
            b = m;
    }

    return b + 2;
}
void e_BSR(BigFloat &P,BigFloat &Q,uint32_t a,uint32_t b,int tds = 1){
    //  Binary Splitting recursion for exp(1).

    if (b - a == 1){
        P = BigFloat(1);
        Q = BigFloat(b);
        return;
    }

    uint32_t m = (a + b) / 2;

    BigFloat P0,Q0,P1,Q1;

    if (b - a < 1000 || tds < 2){
        //  No more threads.
        e_BSR(P0,Q0,a,m);
        e_BSR(P1,Q1,m,b);
    }else{
        //  Run sub-recursions in parallel.
        int tds0 = tds / 2;
        int tds1 = tds - tds0;
#pragma omp parallel num_threads(2)
        {
            int tid = omp_get_thread_num();
            if (tid == 0){
                e_BSR(P0,Q0,a,m,tds0);
            }
            if (tid != 0 || omp_get_num_threads() < 2){
                e_BSR(P1,Q1,m,b,tds1);
            }
        }
    }

    P = P0.mul(Q1,0,tds).add(P1);
    Q = Q0.mul(Q1,0,tds);
}
void e(size_t digits,int tds){
    //  The leading 2 doesn't count.
    digits++;

    size_t p = (digits + 8) / 9;
    size_t terms = e_terms(p);

    //  Limit Exceeded
    if ((uint32_t)terms != terms)
        throw "Limit Exceeded";

    cout << "Computing e..." << endl;
    cout << "Algorithm: Taylor Series of exp(1)" << endl << endl;

    double time0 = wall_clock();

    cout << "Summing Series... " << terms << " terms" << endl;
    BigFloat P,Q;
    e_BSR(P,Q,0,(uint32_t)terms,tds);
    double time1 = wall_clock();
    cout << "Time: " << time1 - time0 << endl;
    
    cout << "Division... " << endl;
    P = P.div(Q,p,tds).add(BigFloat(1),p);
    double time2 = wall_clock();
    cout << "Time: " << time2 - time1 << endl;

    cout << "Total Time = " << time2 - time0 << endl << endl;

    dump_to_file("e.txt",P.to_string(digits));
}
}

int main(){
    omp_set_nested(1);
    int threads = omp_get_max_threads();

    size_t digits = 1000000;

    //  Figure out how large to make the table:
    //  Determine minimum FFT size.
    size_t p = 2*digits / 9 + 10;
    int k = 0;
    size_t length = 1;
    while (length < 3*p){
        length <<= 1;
        k++;
    }

    Bluecrunch::fft_ensure_table(k);
    Bluecrunch::e (digits,threads);
}
