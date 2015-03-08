#define FFT_THRESHOLD_K     20

typedef struct my_complex{
  double r;
  double i;
} my_complex_t;

extern my_complex_t* twiddle_table[32];
extern int twiddle_table_size;

void fft_ensure_table(int k);
void fft_forward(__m128d *T,int k,int tds);
void fft_inverse(__m128d *T,int k,int tds);
void fft_forward_uncached(__m128d *T,int k,int tds);
void fft_inverse_uncached(__m128d *T,int k,int tds);
void fft_pointwise(__m128d *T,__m128d *A,int k);
void int_to_fft(__m128d *T,int k,const uint32_t *A,size_t AL,int digits_per_point);
void fft_to_int(__m128d *T,int k,uint32_t *A,size_t AL,int digits_per_point);