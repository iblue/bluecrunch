#include <math.h>
#include "fft.h"

void fftnew_forward(complex double* T, size_t len) {
  if(len == 1) {
    return;
  }

  //F_2 = [ 1  1 ]
  //      [ 1 -1 ]
  if(len == 2) {
    complex double a = T[0];
    complex double b = T[1];
    T[0] = a+b;
    T[1] = a-b;
    return;
  }

  //       [ 1   1   1 ]
  // F_3 = [ 1 w^1 w^2 ] with w = exp(2*pi*i/3)
  //       [ 1 w^2 w^1 ]
  if(len == 3) {
    complex double omega_1_3 = -0.5 + 0.8660254037844386467637231707529361834714026269051903140*I;
    complex double omega_2_3 = -0.5 - 0.8660254037844386467637231707529361834714026269051903140*I;

    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];

    T[0] = a + b + c;
    T[1] = a + omega_1_3 * b + omega_2_3 * c;
    T[2] = a + omega_2_3 * b + omega_1_3 * c;
    return;
  }

  //       [ 1   1   1   1 ]
  // F_4 = [ 1   i  -1  -i ]
  //       [ 1  -1   1  -1 ]
  //       [ 1  -i  -1   i ]
  if(len == 4) {
    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];
    complex double d = T[3];

    T[0] = a + b + c + d;
    T[1] = a + I*b -c -I*d;
    T[2] = a - b + c - d;
    T[3] = a - I*b -c +I*d;
    return;
  }

  // F_8 = ...
  /*
  if(len == 8) {
    double invsqrt_2 = 0.707106781186547524400844362104849039284835937688474036588;
    complex double W[] = {
      1,
      invsqrt_2 + invsqrt_2*I,
      I,
      -invsqrt_2 + invsqrt_2*I,
      -1,
      -invsqrt_2 - invsqrt_2*I,
      -I,
      invsqrt_2 - invsqrt_2*I
    };

    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];
    complex double d = T[3];
    complex double e = T[4];
    complex double f = T[5];
    complex double g = T[6];
    complex double h = T[7];

    T[0] = a +      b +      c +      d +      e +      f +      g +      h;
    T[1] = a + W[1]*b + W[2]*c + W[3]*d + W[4]*e + W[5]*f + W[6]*g + W[7]*h;
    T[2] = a + W[2]*b + W[4]*c + W[6]*d +      e + W[2]*f + W[4]*g + W[6]*h;
    T[3] = a + W[3]*b + W[6]*c + W[1]*d + W[4]*e + W[7]*f + W[2]*g + W[5]*h;
    T[4] = a + W[4]*b +      c + W[4]*d +      e + W[4]*f +      g + W[4]*h;
    T[5] = a + W[5]*b + W[2]*c + W[7]*d + W[4]*e + W[1]*f + W[6]*g + W[3]*h;
    T[6] = a + W[6]*b + W[4]*c + W[2]*d +      e + W[6]*f + W[4]*g + W[2]*h;
    T[7] = a + W[7]*b + W[6]*c + W[5]*d + W[4]*e + W[3]*f + W[2]*g + W[1]*h;
    return;
  }
  */

  {
    double omega = 2*3.141592653589793238462643383279502884197169399375105820974/len;

    size_t quarter_len = len/4;

    // Radix-4 reduction
    //
    for(size_t k=0;k<quarter_len;k++) {
      // FIXME: Use table
      complex double tf1 = cos(  omega*k) + I*sin(  omega*k);
      complex double tf2 = cos(2*omega*k) + I*sin(2*omega*k);
      complex double tf3 = cos(3*omega*k) + I*sin(3*omega*k);

      complex double a = T[k];
      complex double b = T[k+quarter_len];
      complex double c = T[k+2*quarter_len];
      complex double d = T[k+3*quarter_len];

      T[k]               =     (a +   b + c +   d);
      T[k+  quarter_len] = tf1*(a + I*b - c - I*d);
      T[k+2*quarter_len] = tf2*(a -   b + c -   d);
      T[k+3*quarter_len] = tf3*(a - I*b - c + I*d);
    }

    // Run FFT on the quarters
    // FIXME: Parallel
    fftnew_forward(T,               quarter_len);
    fftnew_forward(T+  quarter_len, quarter_len);
    fftnew_forward(T+2*quarter_len, quarter_len);
    fftnew_forward(T+3*quarter_len, quarter_len);
  }
}

static inline void _fftnew_forward_recursion(complex double* T, size_t len) {
  return;
}

void fftnew_inverse(complex double* T, size_t len) {
  if(len == 1) {
    return;
  }

  //F_2^H = [ 1  1 ]
  //      [ 1 -1 ]
  if(len == 2) {
    complex double a = T[0];
    complex double b = T[1];
    T[0] = (a+b)/2.0;
    T[1] = (a-b)/2.0;
    return;
  }

  //         [ 1   1   1 ]
  // F_3^H = [ 1 w^1 w^2 ] with w = exp(2*pi*(-i)/3)
  //         [ 1 w^2 w^1 ]
  if(len == 3) {
    complex double omega_1_3 = -0.5 - 0.8660254037844386467637231707529361834714026269051903140*I;
    complex double omega_2_3 = -0.5 + 0.8660254037844386467637231707529361834714026269051903140*I;

    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];

    T[0] = (a + b + c)/3.0;
    T[1] = (a + omega_1_3 * b + omega_2_3 * c)/3.0;
    T[2] = (a + omega_2_3 * b + omega_1_3 * c)/3.0;
  }

  //         [ 1   1   1   1 ]
  // F_4^H = [ 1  -i  -1   i ]
  //         [ 1  -1   1  -1 ]
  //         [ 1   i  -1  -i ]
  if(len == 4) {
    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];
    complex double d = T[3];

    T[0] = (a + b + c + d)/4.0;
    T[1] = (a - I*b -c +I*d)/4.0;
    T[2] = (a - b + c - d)/4.0;
    T[3] = (a + I*b -c -I*d)/4.0;
    return;
  }

  // F_8^H = ...
  /*
  if(len == 8) {
    double invsqrt_2 = 0.707106781186547524400844362104849039284835937688474036588;
    complex double W[] = {
      1,
      invsqrt_2 - invsqrt_2*I,
      I,
      -invsqrt_2 - invsqrt_2*I,
      -1,
      -invsqrt_2 + invsqrt_2*I,
      -I,
      invsqrt_2 + invsqrt_2*I
    };

    complex double a = T[0];
    complex double b = T[1];
    complex double c = T[2];
    complex double d = T[3];
    complex double e = T[4];
    complex double f = T[5];
    complex double g = T[6];
    complex double h = T[7];

    T[0] = (a +      b +      c +      d +      e +      f +      g +      h)/8.0;
    T[1] = (a + W[1]*b + W[2]*c + W[3]*d + W[4]*e + W[5]*f + W[6]*g + W[7]*h)/8.0;
    T[2] = (a + W[2]*b + W[4]*c + W[6]*d +      e + W[2]*f + W[4]*g + W[6]*h)/8.0;
    T[3] = (a + W[3]*b + W[6]*c + W[1]*d + W[4]*e + W[7]*f + W[2]*g + W[5]*h)/8.0;
    T[4] = (a + W[4]*b +      c + W[4]*d +      e + W[4]*f +      g + W[4]*h)/8.0;
    T[5] = (a + W[5]*b + W[2]*c + W[7]*d + W[4]*e + W[1]*f + W[6]*g + W[3]*h)/8.0;
    T[6] = (a + W[6]*b + W[4]*c + W[2]*d +      e + W[6]*f + W[4]*g + W[2]*h)/8.0;
    T[7] = (a + W[7]*b + W[6]*c + W[5]*d + W[4]*e + W[3]*f + W[2]*g + W[1]*h)/8.0;
  }
  */

  {
    size_t quarter_len = len/4;

    // Run FFT on the quarters
    fftnew_inverse(T,               quarter_len);
    fftnew_inverse(T+  quarter_len, quarter_len);
    fftnew_inverse(T+2*quarter_len, quarter_len);
    fftnew_inverse(T+3*quarter_len, quarter_len);

    double omega = 2*3.141592653589793238462643383279502884197169399375105820974/len;


    // Radix-4 reduction
    //
    for(size_t k=0;k<quarter_len;k++) {
      // FIXME: Use table
      complex double tf1 = cos(  omega*k) + I*sin(  omega*k);
      complex double tf2 = cos(2*omega*k) + I*sin(2*omega*k);
      complex double tf3 = cos(3*omega*k) + I*sin(3*omega*k);

      complex double a = T[k];
      complex double b = T[k+quarter_len]   * tf1;
      complex double c = T[k+2*quarter_len] * tf2;
      complex double d = T[k+3*quarter_len] * tf3;

      T[k]               = (a +   b + c +   d)/4.0;
      T[k+  quarter_len] = (a - I*b - c + I*d)/4.0;
      T[k+2*quarter_len] = (a -   b + c -   d)/4.0;
      T[k+3*quarter_len] = (a + I*b - c - I*d)/4.0;
    }
  }
}
