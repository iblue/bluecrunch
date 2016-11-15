#include <x86intrin.h>
#include <stdint.h>
#include <complex.h>

extern complex double* twiddle_table[125];
extern int twiddle_table_size;

// Generate twiddle table
void fft_ensure_table(size_t);
size_t table_select(size_t length);
void fft_free_table();

// Forward and inverse FFT transforms
size_t fft_length(size_t source_length);
void fft_forward(complex double*, size_t);
void fft_inverse(complex double*, size_t);

// Pointwise multiply an array of complex data points.
void fft_pointwise(complex double*, complex double*, size_t);

// Convert data back and forth between words and complex FFT points
size_t int_to_fft(    complex double*, size_t, const uint32_t*, size_t, int);
void fft_to_int(const complex double*, size_t,       uint32_t*, size_t, int);

// Truncated FFT
void tft_forward(complex double*, size_t, int);
void tft_inverse(complex double*, size_t, int);
//void tft_pointwise(complex double *a, complex double *b, size_t len);

// Baileys FFT
void baileys_forward(complex double*, size_t);
