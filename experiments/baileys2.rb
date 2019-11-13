# Baileys Algorithmn Prototype

require "byebug"

# Gopnik bit reverse (extra fast)
def bitreverse(int, size)
  len = Math.log2(size)
  int.to_s(2).rjust(len, '0').reverse.to_i(2)
end

# Irish cream FFT
def baileys(values)
  b = values.size
  a = 1
  while a < b
    a*=2
    b/=2
  end

  # Example:
  # a = 8
  # b = 4
  #
  # Before:
  #
  #
  #     b=4 cols
  #
  #  0  1  2  3
  #  4  5  6  7
  #  8  9 10 11
  # 12 13 14 15
  # 16 17 18 19 a=8 rows
  # 20 21 22 23
  # 24 25 26 27
  # 28 29 30 31

  # Do b a-point ffts on the cols (in-place, strided)
  puts "Doing #{b}x #{a}-point FFT in the rows"
  b.times do |i|
    values = stridedfft(values, b, i)
    puts values.to_s
  end

  # Afterwards (columns are now bit reverse)
  # 0 1 2 3 4 5 6 7 -> 0 4 2 6 1 5 3 7
  #
  #    b=4 cols
  #
  #  0  1  2  3
  # 16 17 18 19
  #  8  9 10 11
  # 24 25 26 27 a=8 rows
  #  4  5  6  7
  # 20 21 22 23
  # 12 13 14 15
  # 28 29 30 31

  puts values

  # Multiply with twiddles in-place (twiddles need to be bit-reversed in the cols)
  # maybe merged with the strided fft.
  #
  omega = 2*Math::PI/values.size
  a.times do |i|
    # Get correct row index
    istar = bitreverse(i ,a)

    b.times do |j|
      # a = 8 (index i)
      # b = 4 (index j)

      k = istar*j
      angle = omega*k
      real = Math.cos(angle)
      imag = Math.sin(angle)
      t = Complex(real, imag)

      puts "Matrix #{i},#{j} (#{i*b+j}) (which is really #{istar},#{j} (#{istar*b+j}), twiddle: #{t}"

      values[i*b+j] *= t
    end
  end

  puts values

  # Now do the FFT in the rows (will be very fast :)
  puts "Doing #{a}x #{b}-point FFT in the rows"
  a.times do |i|
    values[i*b...(i+1)*b] = fft(values[i*b...(i+1)*b])
  end

  # Ordering afterwards:
  # (which is basically the same as the normal fft, but transposed)
  #
  # 0 1 2 3 -> 0 2 1 3
  #
  #    b=4 cols
  #
  #  0  2  1  3
  # 16 18 17 19
  #  8 10  9 11
  # 24 26 25 27 a=8 rows
  #  4  6  5  7
  # 20 22 21 23
  # 12 14 13 15
  # 28 30 29 31
  #

  values
end

def ibaileys(values)
  b = values.size
  a = 1
  while a < b
    a*=2
    b/=2
  end
  a, b = [b, a] # works equally well


  puts "Doing #{a}x #{b}-point iFFT in the rows"
  a.times do |i|
    values[i*b...(i+1)*b] = ifft(values[i*b...(i+1)*b])
  end

  # Multiply with twiddles in-place (twiddles need to be bit-reversed in the cols)
  # maybe merged with the strided ifft?
  #
  omega = -2*Math::PI/values.size
  a.times do |i|
    # Get correct row index
    istar = bitreverse(i ,a)

    b.times do |j|
      # a = 8 (index i)
      # b = 4 (index j)

      k = istar*j
      angle = omega*k
      real = Math.cos(angle)
      imag = Math.sin(angle)
      t = Complex(real, imag)

      puts "Matrix #{i},#{j} (#{i*b+j}) (which is really #{istar},#{j} (#{istar*b+j}), twiddle: #{t}"

      values[i*b+j] *= t
    end
  end

  # Do b a-point ffts on the cols (in-place, strided)
  puts "Doing #{b}x #{a}-point FFT in the rows"
  b.times do |i|
    values = stridedifft(values, b, i)
    puts values.to_s
  end

  values
end

def fft(values)
  n = values.length
  #values = values.dup

  if n == 1
    return values
  end

  omega = 2*Math::PI/n

  (0..n/2-1).each do |k|
    # Calc twidles
    angle = omega*k
    r = Math.cos(angle)
    i = Math.sin(angle)

    t = Complex(r, i)

    a = values[k]
    b = values[k+n/2]

    values[k] = (a+b).to_c
    values[k+n/2] = (a-b)*t
  end

  values[0..n/2-1] = fft(values[0..n/2-1])
  values[n/2..-1]  = fft(values[n/2..-1])

  values
end

def stridedfft(values, stride, shift)
  n = values.length/stride

  puts "  strided fft: #{values.length}, #{stride}, #{shift}"


  if n == 1
    return values
  end

  omega = 2*Math::PI/n

  (0..n/2-1).each do |k|
    # Calc twidles
    angle = omega*k
    r = Math.cos(angle)
    i = Math.sin(angle)

    t = Complex(r, i)

    a = values[k*stride+shift]
    b = values[(k+n/2)*stride+shift]

    values[k*stride+shift] = (a+b).to_c
    values[(k+n/2)*stride+shift] = (a-b)*t
  end

  values[0..(n*stride/2-1)] = stridedfft(values[0..(n*stride/2-1)], stride, shift)
  values[(n*stride/2)..-1]  = stridedfft(values[(n*stride/2)..-1], stride, shift)

  values
end


def stridedifft(values, stride, shift)
  n = values.length/stride

  puts "  strided ifft: #{values.length}, #{stride}, #{shift}"


  if n == 1
    return values
  end

  values[0..(n*stride/2-1)] = stridedifft(values[0..(n*stride/2-1)], stride, shift)
  values[(n*stride/2)..-1]  = stridedifft(values[(n*stride/2)..-1], stride, shift)

  omega = -2*Math::PI/n

  (0..n/2-1).each do |k|
    # Calc twidles
    angle = omega*k
    r = Math.cos(angle)
    i = Math.sin(angle)

    t = Complex(r, i)

    a = values[k*stride+shift]
    b = values[(k+n/2)*stride+shift]*t

    values[k*stride+shift] = a+b
    values[(k+n/2)*stride+shift] = a-b
  end

  values
end

def ifft(values)
  n = values.length

  if n == 1
    return values
  end

  values[0..n/2-1] = ifft(values[0..n/2-1])
  values[n/2..-1]  = ifft(values[n/2..-1])

  omega = -2*Math::PI/n

  (0..n/2-1).each do |k|
    # Calc twidles
    angle = omega*k
    r = Math.cos(angle)
    i = Math.sin(angle)

    t = Complex(r, i)

    a = values[k]
    b = values[k+n/2]*t

    values[k] = a+b
    values[k+n/2] = a-b
  end

  values
end

def int_to_fft(int)
  arr = int.to_s.split(//).map(&:to_i).reverse
  i = 1
  i*=2 while i<arr.size
  i*=2
  ret = Array.new(i,0)
  ret[0..arr.size-1] = arr
  ret
end

def fft_to_int(arr)
  # scale and convert
  arr = arr.map{|x| (x*(1.0/arr.size)).real.round.to_i}

  # carry (gopnik variant)
  c = 1
  ret = arr.map do |v|
    r = v*c
    c*=10
    r
  end
  ret.sum
end

def fft_pointwise(a,b)
  a.zip(b).map{|x,y| x*y}
end

def mul(a,b)
  a = int_to_fft(a)
  b = int_to_fft(b)
  a = baileys(a)
  b = baileys(b)
  c = fft_pointwise(a,b)
  c = ibaileys(c)
  c = fft_to_int(c)
end

if false
  a = 1234541623518273
  b = 1234557781293722
  c = mul(a,b)
end

if false
  a = [1, 5, 3, 6, 17, 4, 18, 22]
  stride = 4
  b = a.dup.inject([]) { |x,y| x + [0] + [y] + [0]*(stride-2)}

  puts a.to_s
  puts b.to_s

  puts "FFT:"
  a = fft(a)
  puts a
  puts "Strided:"
  b = stridedfft(b, stride, 1)
  puts b
  exit
end

if true
  a = fft = [
  775996.000000+0.000000i, # 0
  290480.000000+0.000000i, # 1
  915224.000000+0.000000i, # 2
  558577.000000+0.000000i, # 3
  106561.000000+0.000000i, # 4
  42574.000000+0.000000i, # 5
  726212.000000+0.000000i, # 6
  187723.000000+0.000000i, # 7
  113371.000000+0.000000i, # 8
  56212.000000+0.000000i, # 9
  572015.000000+0.000000i, # 10
  660642.000000+0.000000i, # 11
  825370.000000+0.000000i, # 12
  576169.000000+0.000000i, # 13
  692923.000000+0.000000i, # 14
  885385.000000+0.000000i, # 15
  872504.000000+0.000000i, # 16
  759674.000000+0.000000i, # 17
  45103.000000+0.000000i, # 18
  689804.000000+0.000000i, # 19
  85124.000000+0.000000i, # 20
  999203.000000+0.000000i, # 21
  126354.000000+0.000000i, # 22
  279173.000000+0.000000i, # 23
  249720.000000+0.000000i, # 24
  753803.000000+0.000000i, # 25
  1493.000000+0.000000i, # 26
  177949.000000+0.000000i, # 27
  390365.000000+0.000000i, # 28
  393232.000000+0.000000i, # 29
  876465.000000+0.000000i, # 30
  955656.000000+0.000000i, # 31
  920222.000000+0.000000i, # 32
  1033741.000000+0.000000i, # 33
  741732.000000+0.000000i, # 34
  327372.000000+0.000000i, # 35
  162649.000000+0.000000i, # 36
  973912.000000+0.000000i, # 37
  300865.000000+0.000000i, # 38
  579844.000000+0.000000i, # 39
  941890.000000+0.000000i, # 40
  653221.000000+0.000000i, # 41
  737896.000000+0.000000i, # 42
  850977.000000+0.000000i, # 43
  616796.000000+0.000000i, # 44
  261488.000000+0.000000i, # 45
  1036398.000000+0.000000i, # 46
  930821.000000+0.000000i, # 47
  784064.000000+0.000000i, # 48
  955342.000000+0.000000i, # 49
  1032821.000000+0.000000i, # 50
  57939.000000+0.000000i, # 51
  560737.000000+0.000000i, # 52
  608138.000000+0.000000i, # 53
  297266.000000+0.000000i, # 54
  58129.000000+0.000000i, # 55
  663976.000000+0.000000i, # 56
  717459.000000+0.000000i, # 57
  252517.000000+0.000000i, # 58
  853719.000000+0.000000i, # 59
  153208.000000+0.000000i, # 60
  516130.000000+0.000000i, # 61
  50259.000000+0.000000i, # 62
  186645.000000+0.000000i, # 63
  303975.000000+0.000000i, # 64
  597650.000000+0.000000i, # 65
  995580.000000+0.000000i, # 66
  653022.000000+0.000000i, # 67
  865239.000000+0.000000i, # 68
  792062.000000+0.000000i, # 69
  807097.000000+0.000000i, # 70
  491531.000000+0.000000i, # 71
  299103.000000+0.000000i, # 72
  551239.000000+0.000000i, # 73
  371809.000000+0.000000i, # 74
  218552.000000+0.000000i, # 75
  1012834.000000+0.000000i, # 76
  854134.000000+0.000000i, # 77
  922650.000000+0.000000i, # 78
  14186.000000+0.000000i, # 79
  363702.000000+0.000000i, # 80
  694458.000000+0.000000i, # 81
  733849.000000+0.000000i, # 82
  1024241.000000+0.000000i, # 83
  736113.000000+0.000000i, # 84
  73447.000000+0.000000i, # 85
  831046.000000+0.000000i, # 86
  914741.000000+0.000000i, # 87
  668232.000000+0.000000i, # 88
  517310.000000+0.000000i, # 89
  922462.000000+0.000000i, # 90
  331172.000000+0.000000i, # 91
  137106.000000+0.000000i, # 92
  383005.000000+0.000000i, # 93
  952610.000000+0.000000i, # 94
  383219.000000+0.000000i, # 95
  0.000000+0.000000i, # 96
  0.000000+0.000000i, # 97
  0.000000+0.000000i, # 98
  0.000000+0.000000i, # 99
  0.000000+0.000000i, # 100
  0.000000+0.000000i, # 101
  0.000000+0.000000i, # 102
  0.000000+0.000000i, # 103
  0.000000+0.000000i, # 104
  0.000000+0.000000i, # 105
  0.000000+0.000000i, # 106
  0.000000+0.000000i, # 107
  0.000000+0.000000i, # 108
  0.000000+0.000000i, # 109
  0.000000+0.000000i, # 110
  0.000000+0.000000i, # 111
  0.000000+0.000000i, # 112
  0.000000+0.000000i, # 113
  0.000000+0.000000i, # 114
  0.000000+0.000000i, # 115
  0.000000+0.000000i, # 116
  0.000000+0.000000i, # 117
  0.000000+0.000000i, # 118
  0.000000+0.000000i, # 119
  0.000000+0.000000i, # 120
  0.000000+0.000000i, # 121
  0.000000+0.000000i, # 122
  0.000000+0.000000i, # 123
  0.000000+0.000000i, # 124
  0.000000+0.000000i, # 125
  0.000000+0.000000i, # 126
  0.000000+0.000000i, # 127
  0.000000+0.000000i, # 128
  0.000000+0.000000i, # 129
  0.000000+0.000000i, # 130
  0.000000+0.000000i, # 131
  0.000000+0.000000i, # 132
  0.000000+0.000000i, # 133
  0.000000+0.000000i, # 134
  0.000000+0.000000i, # 135
  0.000000+0.000000i, # 136
  0.000000+0.000000i, # 137
  0.000000+0.000000i, # 138
  0.000000+0.000000i, # 139
  0.000000+0.000000i, # 140
  0.000000+0.000000i, # 141
  0.000000+0.000000i, # 142
  0.000000+0.000000i, # 143
  0.000000+0.000000i, # 144
  0.000000+0.000000i, # 145
  0.000000+0.000000i, # 146
  0.000000+0.000000i, # 147
  0.000000+0.000000i, # 148
  0.000000+0.000000i, # 149
  0.000000+0.000000i, # 150
  0.000000+0.000000i, # 151
  0.000000+0.000000i, # 152
  0.000000+0.000000i, # 153
  0.000000+0.000000i, # 154
  0.000000+0.000000i, # 155
  0.000000+0.000000i, # 156
  0.000000+0.000000i, # 157
  0.000000+0.000000i, # 158
  0.000000+0.000000i, # 159
  0.000000+0.000000i, # 160
  0.000000+0.000000i, # 161
  0.000000+0.000000i, # 162
  0.000000+0.000000i, # 163
  0.000000+0.000000i, # 164
  0.000000+0.000000i, # 165
  0.000000+0.000000i, # 166
  0.000000+0.000000i, # 167
  0.000000+0.000000i, # 168
  0.000000+0.000000i, # 169
  0.000000+0.000000i, # 170
  0.000000+0.000000i, # 171
  0.000000+0.000000i, # 172
  0.000000+0.000000i, # 173
  0.000000+0.000000i, # 174
  0.000000+0.000000i, # 175
  0.000000+0.000000i, # 176
  0.000000+0.000000i, # 177
  0.000000+0.000000i, # 178
  0.000000+0.000000i, # 179
  0.000000+0.000000i, # 180
  0.000000+0.000000i, # 181
  0.000000+0.000000i, # 182
  0.000000+0.000000i, # 183
  0.000000+0.000000i, # 184
  0.000000+0.000000i, # 185
  0.000000+0.000000i, # 186
  0.000000+0.000000i, # 187
  0.000000+0.000000i, # 188
  0.000000+0.000000i, # 189
  0.000000+0.000000i, # 190
  0.000000+0.000000i, # 191
  0.000000+0.000000i, # 192
  0.000000+0.000000i, # 193
  0.000000+0.000000i, # 194
  0.000000+0.000000i, # 195
  0.000000+0.000000i, # 196
  0.000000+0.000000i, # 197
  0.000000+0.000000i, # 198
  0.000000+0.000000i, # 199
  0.000000+0.000000i, # 200
  0.000000+0.000000i, # 201
  0.000000+0.000000i, # 202
  0.000000+0.000000i, # 203
  0.000000+0.000000i, # 204
  0.000000+0.000000i, # 205
  0.000000+0.000000i, # 206
  0.000000+0.000000i, # 207
  0.000000+0.000000i, # 208
  0.000000+0.000000i, # 209
  0.000000+0.000000i, # 210
  0.000000+0.000000i, # 211
  0.000000+0.000000i, # 212
  0.000000+0.000000i, # 213
  0.000000+0.000000i, # 214
  0.000000+0.000000i, # 215
  0.000000+0.000000i, # 216
  0.000000+0.000000i, # 217
  0.000000+0.000000i, # 218
  0.000000+0.000000i, # 219
  0.000000+0.000000i, # 220
  0.000000+0.000000i, # 221
  0.000000+0.000000i, # 222
  0.000000+0.000000i, # 223
  0.000000+0.000000i, # 224
  0.000000+0.000000i, # 225
  0.000000+0.000000i, # 226
  0.000000+0.000000i, # 227
  0.000000+0.000000i, # 228
  0.000000+0.000000i, # 229
  0.000000+0.000000i, # 230
  0.000000+0.000000i, # 231
  0.000000+0.000000i, # 232
  0.000000+0.000000i, # 233
  0.000000+0.000000i, # 234
  0.000000+0.000000i, # 235
  0.000000+0.000000i, # 236
  0.000000+0.000000i, # 237
  0.000000+0.000000i, # 238
  0.000000+0.000000i, # 239
  0.000000+0.000000i, # 240
  0.000000+0.000000i, # 241
  0.000000+0.000000i, # 242
  0.000000+0.000000i, # 243
  0.000000+0.000000i, # 244
  0.000000+0.000000i, # 245
  0.000000+0.000000i, # 246
  0.000000+0.000000i, # 247
  0.000000+0.000000i, # 248
  0.000000+0.000000i, # 249
  0.000000+0.000000i, # 250
  0.000000+0.000000i, # 251
  0.000000+0.000000i, # 252
  0.000000+0.000000i, # 253
  0.000000+0.000000i, # 254
  0.000000+0.000000i, # 255
]

  original = a.dup
  b = a.dup

  puts "FFT:"
  a = fft(a)
  puts a
  puts " === Baileys === "
  b = baileys(b)
  puts b
  puts "inverse FFT:"
  c = ifft(b)
  c = c.map{|x| (x*(1.0/c.size)).real.round.to_i}
  puts c.to_s
  if c == original
    puts "GREAT SUCCESS"
  else
    puts "FAIL. #{c.to_s}"
  end
end

if a*b == c
  puts "GREAT SUCCESS"
else
  puts "FAIL #{a}*#{b} = #{c} (expected: #{a*b})"
end
