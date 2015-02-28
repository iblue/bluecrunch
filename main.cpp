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
#include "bigfloat.h"

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
