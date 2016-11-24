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

  // Korn-Lambiotte
  if(len <= 256) {
    // Recursive splitting
  }
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

  // Korn-Lambiotte
  if(len <= 256) {
    // Recursive splitting
  }
}
