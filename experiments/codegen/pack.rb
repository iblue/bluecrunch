require "byebug"

def bitmask(bits)
  "0x#{((2**bits)-1).to_s(16)}"
end

def wordpos(bits)
  bits/32
end

def missing_bits(bitpos)
  32 - bitpos%32
end

def gencode(bits_per_point)
  # First find out how many words we need to get back to aligning.
  total_bits = bits_per_point.lcm(32)
  align = total_bits/32

  puts "// need alignment #{align}"

  puts "static inline void int_to_fft#{bits_per_point}(compley double *V, const uint32_t *A, size_t AL) {"
  puts "  __m128d* T = (__m128*)V;"
  puts "  complex double *origV = V;"
  puts ""
  puts "  for(size_t c = 0; c < AL/#{align}+#{align-1}; c++) {"

  bitpos  = 0
  while true do
    # We crossed a word boundary
    wordpos = wordpos(bitpos)
    if(wordpos(bitpos-bits_per_point) != wordpos)
      puts "    if(#{align}*c+#{wordpos} >= AL) {"
      if missing_bits(bitpos-bits_per_point) && wordpos != 0
        puts "      *T++ = _mm_set_sd(word#{wordpos-1} & #{bitmask(missing_bits(bitpos-bits_per_point))});"
      end
      puts "      break;"
      puts "    }"
      puts ""
      puts "    uint32_t word#{wordpos} = A[#{align}*c+#{wordpos}];"
      puts ""
    end

    # Get mb bits from the current word get the rest from next
    mb = [missing_bits(bitpos), bits_per_point].min
    if mb == bits_per_point
      puts "    *T++ = _mm_set_sd(word#{wordpos} & #{bitmask(mb)});"
      puts "    word#{wordpos} >>= #{mb};"
    else
      rb = bits_per_point - mb
      puts "    *T++ = _mm_set_sd((word#{wordpos-1} & #{bitmask(mb)}) | (word#{wordpos} & #{bitmask(rb)});"
      puts "    word#{wordpos} >>= #{mb};"
    end

    bitpos += bits_per_point

    return if total_bits == bitpos
  end
  puts "  }"
  puts "}"
end

#(8..22).each do |i|
#  puts "// Code for #{i} bits per point"
#  gencode(i)
#end

gencode(12)
