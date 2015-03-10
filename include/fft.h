#define FFT_THRESHOLD_K     20

typedef struct my_complex {
  double r;
  double i;
} my_complex_t;

extern my_complex_t* twiddle_table[32];
extern int twiddle_table_size;

// Generate twiddle table
void fft_ensure_table(int);

// Forward and inverse FFT transforms
void fft_forward(__m128d*, int, int);
void fft_inverse(__m128d*, int, int);
void fft_forward_uncached(__m128d*, int, int);
void fft_inverse_uncached(__m128d*, int, int);

// Pointwise multiply an array of complex data points.
void fft_pointwise(__m128d*, __m128d *, int);

// Convert data back and forth between words and complex FFT points
void int_to_fft(      __m128d*, int, const uint32_t*, size_t, int);
void fft_to_int(const __m128d*, int,       uint32_t*, size_t, int);
