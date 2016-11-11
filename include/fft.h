#include <x86intrin.h>
#include <stdint.h>
#include <complex.h>

// When do we start to decompose FFTs for parallization
#define FFT_THRESHOLD_K 10

extern complex double* twiddle_table[32];
extern int twiddle_table_size;

// Generate twiddle table
void fft_ensure_table(int);

// Forward and inverse FFT transforms
void fft_forward(complex double*, int);
void fft_inverse(complex double*, int);
void fft_forward_uncached(complex double*, int);
void fft_inverse_uncached(complex double*, int);

// Pointwise multiply an array of complex data points.
void fft_pointwise(complex double*, complex double*, int);

// Convert data back and forth between words and complex FFT points
size_t int_to_fft(    complex double*, int, const uint32_t*, size_t, int);
void fft_to_int(const complex double*, int,       uint32_t*, size_t, int);

// Truncated FFT
void tft_forward(complex double*, size_t, int);
void tft_inverse(complex double*, size_t, int);
//void tft_pointwise(complex double *a, complex double *b, size_t len);

// Baileys FFT
void baileys_forward(complex double*, int);
