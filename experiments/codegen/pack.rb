require "byebug"

def outcode(code)
  $outfile.write(code + "\n")
end

def bitmask(bits)
  "0x#{((2**bits)-1).to_s(16)}"
end

def init(bpp)
  puts "init #{bpp}"
  code(:init, bpp)
end

def deinit
  puts "ret"
  puts ""
  code(:ret)
end

# Load {bits} bits from {word} word
def get(bits, word, bitpos)
  puts "get #{bits} from #{word} (#{bitpos})"
  code(:get, bits, word, bitpos)
end

# Load word if exists
def ldf(word)
  puts "load_or_fin #{word}"
  code(:ldf, word)
end

def ld(word)
  puts "load #{word}"
  code(:ld, word)
end

def carry(bits)
  puts "carry #{bits}"
  code(:carry, bits)
end

def put
  puts "put"
  code(:put)
end

def conv(pwcnt, pwpos, cwcnt, cwpos)
  if pwcnt == cwcnt
    get(cwpos-pwpos, cwcnt, pwpos)
    put
  else
    bits = 32*cwcnt + cwpos - 32*pwcnt - pwpos

    get(32-pwpos, pwcnt, pwpos)
    if cwpos%bits != 0
      ldf(cwcnt)
      carry(cwpos)
      get(cwpos%bits, cwcnt, cwpos)
    end
    put
  end
end

def gencode(bits_per_point)
  # get 12 bit from word 0 => 0 (20 left)
  # get 12 bit from word 0 => 1 ( 8 left)
  # get  8 bit from word 0 => 2 ( 0 left)
  # get  4 bit from word 1 => 2 (
  # get 12 bit from word 1 => 3
  # get 12 bit from word 1

  bitpos  = 0
  init(bits_per_point)
  ld(0)
  while true do
    pwcnt    = bitpos/32
    pwpos    = bitpos%32
    bitpos += bits_per_point
    cwcnt    = bitpos/32
    cwpos    = bitpos%32

    #puts "#{bitpos}: (#{pwcnt},#{pwpos}) (#{cwcnt},#{cwpos})"
    conv(pwcnt, pwpos, cwcnt, cwpos)

    if cwpos == 0
      break
    end
  end
  deinit
end

def int_to_fft_code(i)
  @ldf = false
  @lg  = 0
  def code(sym, *args)
    case sym
    when :init
      bpp = args[0]
      outcode "// Outgenerated"
      outcode "static inline size_t int_to_fft#{bpp}(complex double *F, const uint32_t *W, size_t WL) {"
      outcode "  uint32_t w=0;"
      outcode "  uint32_t c=0;"
      outcode "  __m128d* T = (__m128d*)F;"
      outcode "  const uint32_t *end = W + WL;"
      outcode "  while(1) {"
    when :ld
      outcode "    w = *W;"
    when :ldf
      outcode "    if(++W >= end) {"
      outcode "      *T++ = _mm_set_sd(c);"
      outcode "      break;"
      outcode "    }"
      outcode "    w = *W;"
      @ldf = true
    when :get
      bits, word, bitpos = args
      if @ldf
        outcode "    c |= (w & #{bitmask(bits)}) << #{@lg};"
        outcode "    w >>= #{bits};"
        @ldf = false
      else
        #outcode "    c <<= #{bits};"
        outcode "    c  |= w & #{bitmask(bits)};"
        outcode "    w >>= #{bits};"
      end
      @lg = bits
    when :put
      outcode "    *T++ = _mm_set_sd(c);"
      outcode "    c = 0;"
    when :ret
      outcode "    if(W++ > end) {"
      outcode "      break;"
      outcode "    }"
      outcode "  }"
      outcode "  return (complex double*)T - F;"
      outcode "}"
      outcode ""
    end
  end
  gencode(i)
end

def fft_to_int_code(i)
  @ldf = false
  def code(sym, *args)
    case sym
    when :init
      bpp = args[0]
      outcode "// Outgenerated"
      outcode "static inline void fft_to_int#{bpp}(const complex double *F, uint32_t *W, size_t WL, double scale) {"
      outcode "  uint64_t c=0;"
      outcode "  uint32_t *end = W + WL;"
      outcode "  double f;"
      outcode "  double i;"
      outcode "  uint32_t w=0;"
      outcode "  while(1) {"
    when :ldf
      outcode "    *W = w;"
      outcode "    if(++W >= end) {"
      outcode "      break;"
      outcode "    }"
      outcode "    w = 0;"
      @ldf = true
    when :get
      bits, word, bitpos = args
      if @ldf
        outcode "    w |= c & #{bitmask(bits)};"
        outcode "    c >>= #{bits};"
        @ldf = false
      else
        outcode "    f  = ((double*)F++)[0] * scale;"
        outcode "    i  = (uint64_t)(f+0.5);"
        outcode "    c += i;"
        outcode "    w |= (c&#{bitmask(bits)}) << #{bitpos};" # Shift into right pos
        outcode "    c >>= #{bits};"
      end
    when :ret
      outcode "    *W = w;"
      outcode "    if(W++ > end) {"
      outcode "      break;"
      outcode "    }"
      outcode "    w = 0;"
      outcode "  }"
      outcode "}"
      outcode ""
    end
  end
  gencode(i)
end


(7..22).each do |i|
  $outfile = File.open("./int_to_fft#{i}.c", "w")
  puts "== Generating Code for #{i}"
  int_to_fft_code(i)
  $outfile = File.open("./fft_to_int#{i}.c", "w")
  fft_to_int_code(i)
end
