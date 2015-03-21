#include <complex.h>
#include <math.h>
#include <x86intrin.h>
#include <stdint.h>
#include "fft.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

/**
 * Implementation of cache friendly tft
 */

static inline int bitlog2(int N) {
  int k = 0;
  while (N >>= 1) {
    k++;
  }
  return k;
}

complex double omega(int i, int N) {
  int k = bitlog2(N);

  complex double* local_table = twiddle_table[k];
  complex double  val         = local_table[i];

  return val;
}

void tft1(complex double *x, int l, int L, int z, int n, int u, int s) {
  if(L == 2) {
    if(n == 2 && z == 2) {
      complex double a = x[u];
      complex double b = x[u+s];
      complex double xi = omega(u%(L*s), L*s);
      x[u]   = a + b;
      x[u+s] = xi*(a - b);
    } else if(n == 2 && z == 1) {
      complex double xi = omega(u%(L*s), L*s);
      x[u+s] = xi*x[u];
    } else if(n == 1 && z == 2) {
      x[u] += x[u+s];
    } else {
    }
    return;
  }

  // recursive case
  int l_1      = l/2;
  int l_2      = l - l/2;
  int L_1      = 1 << l_1;
  int L_2      = 1 << l_2;

  int n_2      = n % L_2; // FIXME: Bitshift?
  int n_1      = n/L_2;
  int n_1prime = (n+(L_2-1))/L_2; // round up by adding L_2-1
  int z_2      = z%L_2;
  int z_1      = z/L_2;

  int z_2prime;
  if(z_1 > 0) {
    z_2prime = L_2;
  } else {
    z_2prime = z_2;
  }


  // col transforms
  for(int k=0;k<z_2;k++) {
    tft1(x, l_1, L_1, z_1+1, n_1prime, u+s*k, s*L_2);
  }
  for(int k=z_2;k<z_2prime;k++) {
    tft1(x, l_1, L_1, z_1, n_1prime, u+s*k, s*L_2);
  }

  // row transforms
  for(int k=0;k<n_1;k++) {
    tft1(x, l_2, L_2, z_2prime, L_2, u+s*k*L_2, s);
  }
  if(n_2 > 0) {
    tft1(x, l_2, L_2, z_2prime, n_2, u+s*n_1*L_2, s);
  }
}

void itft1(complex double *x, int l, int L, int z, int n, int f, int u, int s) {
  if(L == 2) {
    if(n == 2) {
      complex double xi = omega(u%(L*s), L*s);
      xi = creal(xi)-I*cimag(xi);
      complex double a = x[u] + xi*x[u+s];
      complex double b = x[u] - xi*x[u+s];
      x[u]   = a;
      x[u+s] = b;
    }
    if(n == 1 && f == 1 && z == 2) {
      complex double xi = omega(u%(L*s), L*s);
      complex double a = 2*x[u] - x[u+s];
      complex double b = xi*(x[u]-x[u+s]);
      x[u]   = a;
      x[u+s] = b;
    }
    if(n == 1 && f == 1 && z == 1) {
      complex double xi = omega(u%(L*s), L*s);
      complex double a = 2*x[u];
      complex double b = xi*x[u];
      x[u]   = a;
      x[u+s] = b;
    }
    if(n == 1 && f == 0 && z == 2) {
      x[u] = 2*x[u] - x[u+s];
    }
    if(n == 1 && f == 0 && z == 1) {
      x[u] = 2*x[u];
    }
    if(n == 0 && z == 2) {
      x[u] = (x[u] + x[u+s])/2;
    }
    if(n == 0 && z == 1) {
      x[u] = x[u]/2;
    }
    return;
  }

  // recursive case
  int l_1      = l/2;
  int l_2      = l - l/2;
  int L_1      = 1 << l_1;
  int L_2      = 1 << l_2;

  int n_2      = n % L_2; // FIXME: Bitshift?
  int n_1      = n/L_2;
  int z_2      = z%L_2;
  int z_1      = z/L_2;

  int fprime;
  if(n_2 + f > 0) {
    fprime = 1;
  } else {
    fprime = 0;
  }

  int z_2prime;
  if(z_1 > 0) {
    z_2prime = L_2;
  } else {
    z_2prime = z_2;
  }

  int m      = min(n_2, z_2);
  int mprime = max(n_2, z_2);

  // row transforms
  for(int k=0;k<n_1;k++) {
    itft1(x, l_2, L_2, L_2, L_2, 0, u+s*k*L_2, s);
  }

  // rightmost column transforms
  for(int k=n_2;k<mprime;k++) {
    itft1(x, l_1, L_1, z_1+1, n_1, fprime, u+s*k, s*L_2);
  }
  for(int k=mprime;k<z_2prime;k++) {
    itft1(x, l_1, L_1, z_1, n_1, fprime, u+s*k, s*L_2);
  }

  // last row transform
  if(fprime == 1) {
    itft1(x, l_2, L_2, z_2prime, n_2, f, u+s*n_1*L_2, s);
  }

  // leftmost column transforms
  for(int k=0;k<m;k++) {
    itft1(x, l_1, L_1, z_1 + 1, n_1 + 1, 0, u+s*k, s*L_2);
  }
  for(int k=m;k<n_2;k++) {
    itft1(x, l_1, L_1, z_1, n_1 + 1, 0, u+s*k, s*L_2);
  }
}

void tft_forward(complex double *T, int k, size_t in, size_t out) {
  tft1(T, k, 1 << k, in, out, 0, 1);
}

void tft_inverse(complex double *T, int k, size_t in, size_t out) {
  itft1(T, k, 1 << k, in, out, 0, 0, 1);
}
