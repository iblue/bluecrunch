def bitmask(bits)
  "0x#{((2**bits)-1).to_s(16)}"
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
  wordpos = 0
  while true do
    # On word boundary
    puts "    if(#{align}*c+#{wordpos} >= AL) {"
    puts "      break;"
    puts "    }"
    puts ""
    puts "    uint32_t word#{wordpos} = A[#{align}*c+#{wordpos}];"

  end
  puts "  }"
  puts "}"
end

(8..22).each do |i|
  puts "// Code for #{i} bits per point"
  gencode(i)
end
