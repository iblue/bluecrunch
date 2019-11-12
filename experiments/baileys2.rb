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
  a, b = [b, a] # works equally well

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
    b = values[(k+n/2)*stride+shift]

    values[k*stride+shift] = (a+b).to_c
    values[(k+n/2)*stride+shift] = (a-b)*t
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
  puts "FFT input"
  puts a
  puts a.size
  puts "---"
  puts "FFT Output"
  a = fft(a)
  puts a
  puts "---"
  puts "Baileys Output"
  b = baileys(b)
  puts b
  exit
  #b = fft(b)
  #c = fft_pointwise(a,b)
  #c = ifft(c)
  #c = fft_to_int(c)
end

a = 12345
b = 12345
#c = mul(a,b)

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

a = [1, 5, 3, 6, 17, 4, 18, 22, 9, 27, 2, 3, 5, 31, 41, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]
b = a.dup

puts "FFT:"
a = fft(a)
puts a
puts "Baileys:"
b = baileys(b)
puts b
puts "inverse FFT:"
c = ibaileys(b)
c = c.map{|x| (x*(1.0/c.size)).real.round.to_i}
puts c.to_s

exit

if a*b == c
  puts "GREAT SUCCESS"
else
  puts "FAIL #{a}*#{b} = #{c} (expected: #{a*b})"
end
