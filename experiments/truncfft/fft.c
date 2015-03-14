#include <stdio.h>
#include <complex.h>
#include <math.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define max(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a); __typeof__(b) _b = (b); _a < _b ? _a : _b; })

/**
 * Implementation of cache friendly tft
 */

complex double x[16] = {1, 2, 3, 4, 5, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0};

complex double omega(int i, int N) {
  return cexp(2*M_PI*I*1.0/N * i);
}

void tft(int L, complex double xi, int z, int n, int u, int s, int sl) {
  //printf("tft(L=%d, ..., z=%d, n=%d, u=%d, s=%d)\n", L, z, n, u, s);
  if(L == 2) {
    printf("%*scalc n = %d, z = %d\n", sl, "", n, z);
    if(n == 2 && z == 2) {
      complex double a = x[u];
      complex double b = x[u+s];
      x[u]   = a + b;
      x[u+s] = xi*(a - b);
      printf("%*sx[%d], x[%d] = butterfly\n", sl, "", u, u+s);
    } else if(n == 2 && z == 1) {
      x[u+s] = xi*x[u];
      printf("%*sx[%d] = xi*x[%d]\n", sl, "", u+s, u);
    } else if(n == 1 && z == 2) {
      x[u] += x[u+s];
      printf("%*sx[%d] += x[%d]\n", sl, "", u, u+s);
    } else {
      printf("%*sreturn\n", sl, "");
    }
    return;
  }

  int l = log((float)L)/log(2.0); // Cache tuning param FIXME

  // recursive case
  int L_1      = 1 << l/2;
  int L_2      = 1 << (l - l/2);
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

  printf("%*sparams: L = %d, L_1 = %d, L_2 = %d\n", sl, "", L, L_1, L_2);

  // col transforms
  for(int k=0;k<z_2;k++) {
    printf("%*scol offset %d, stride %d (0..z_2)\n", sl, "", u+s*k, s*L_2);
    tft(L_1, omega(k,L)*xi, z_1+1, n_1prime, u+s*k, s*L_2, sl+2);
  }
  for(int k=z_2;k<z_2prime;k++) {
    printf("%*scol offset %d, stride %d (z_2..)\n", sl, "", u+s*k, s*L_2);
    tft(L_1, omega(k,L)*xi, z_1, n_1prime, u+s*k, s*L_2, sl+2);
  }

  // row transforms
  for(int k=0;k<n_1;k++) {
    printf("%*srow offset %d, stride %d (n_1)\n", sl, "", u+s*k*L_2, s);
    tft(L_2, cpow(xi, L_1), z_2prime, L_2, u+s*k*L_2, s, sl+2);
  }
  if(n_2 > 0) {
    printf("%*srow offset %d, stride %d (n_2)\n", sl, "", u+s*n_1*L_2, s);
    tft(L_2, cpow(xi, L_1), z_2prime, n_2, u+s*n_1*L_2, s, sl+2);
  }
}

void itft(int L, complex double xi, int z, int n, int f, int u, int s, int sl) {
  if(L == 2) {
    if(n == 2) {
      complex double a = x[u] + cpow(xi, -1)*x[u+s];
      complex double b = x[u] - cpow(xi, -1)*x[u+s];
      x[u]   = a;
      x[u+s] = b;
    }
    if(n == 1 && f == 1 && z == 2) {
      complex double a = 2*x[u] - x[u+s];
      complex double b = xi*(x[u]-x[u+s]);
      x[u]   = a;
      x[u+s] = b;
    }
    if(n == 1 && f == 1 && z == 1) {
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

  int l = log((float)L)/log(2.0); // Cache tuning param FIXME

  // recursive case
  int L_1      = 1 << l/2;
  int L_2      = 1 << (l - l/2);
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
    itft(L_2, cpow(xi, L_1), L_2, L_2, 0, u+s*k*L_2, s, sl+2);
  }

  // rightmost column transforms
  for(int k=n_2;k<mprime;k++) {
    itft(L_1, omega(k,L)*xi, z_1+1, n_1, fprime, u+s*k, s*L_2, sl+2);
  }
  for(int k=mprime;k<z_2prime;k++) {
    itft(L_1, omega(k,L)*xi, z_1, n_1, fprime, u+s*k, s*L_2, sl+2);
  }

  // last row transform
  if(fprime == 1) {
    itft(L_2, cpow(xi, L_1), z_2prime, n_2, f, u+s*n_1*L_2, s, sl+2);
  }

  // leftmost column transforms
  for(int k=0;k<m;k++) {
    itft(L_1, omega(k,L)*xi, z_1 + 1, n_1 + 1, 0, u+s*k, s*L_2, sl+2);
  }
  for(int k=m;k<n_2;k++) {
    itft(L_1, omega(k,L)*xi, z_1, n_1 + 1, 0, u+s*k, s*L_2, sl+2);
  }
}

void dump() {
  printf("(%f,%f), (%f, %f), (%f, %f), (%f, %f), "
         "(%f,%f), (%f, %f), (%f, %f), (%f, %f)\n",
      creal(x[0]), cimag(x[0]),
      creal(x[1]), cimag(x[1]),
      creal(x[2]), cimag(x[2]),
      creal(x[3]), cimag(x[3]),
      creal(x[4]), cimag(x[4]),
      creal(x[5]), cimag(x[5]),
      creal(x[6]), cimag(x[6]),
      creal(x[7]), cimag(x[7])
  );
}

int main() {
  int N = 8;
  int in = 5;
  int out = 5;
  dump();
  tft(N, 1.0, in, out, 0, 1, 0);
  dump();
  itft(N, 1.0, in, out, 0, 0, 1, 0);
  dump();
}
