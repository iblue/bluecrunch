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
  #  0  1  2  3
  #  4  5  6  7
  #  8  9 10 11
  # 12 13 14 15
  # 16 17 18 19
  # 20 21 22 23
  # 24 25 26 27
  # 28 29 30 31

  # Do a b-point ffts on the rows (in-place)
  puts "Doing #{a}x #{b}-point FFT in the rows"
  a.times do |i|
    values[i*b...(i+1)*b] = fft(values[i*b...(i+1)*b])
  end

  # Afterwards:
  #
  #    b = 4
  #
  #  0  2  1 3
  #  4  6  5 7
  #  8 10  9 11
  # 12 14 13 15
  # 16 18 17 19  a = 8
  # 20 22 21 23
  # 24 26 25 27
  # 28 30 29 31

  puts values

  # Multiply and transpose
  res = Array.new(a*b, 0)
  omega = 2*Math::PI/values.size
  a.times do |i|
    values[i*b...(i+1)*b].each_with_index do |val, j|

      # a = 8 (index i)
      # b = 4 (index j) -> reverse

      # Multiply with twiddles
      jstar = bitreverse(j ,b)
      k = jstar*i # ??? FIXME
      angle = omega*k
      real = Math.cos(angle)
      imag = Math.sin(angle)
      t = Complex(real, imag)
      #t = Complex(1.0, 0.0) # FIXME: Remove

      puts "Index: #{jstar},#{i}: #{val} *===* #{t} (TARGET: #{i+j*a})"

      # TODO Later: movntq
      # FIXME: Do not use jstar. Does not matter, just reorders and makes shit
      # slow (probably). Use j.
      res[i+j*a] = val*t

      puts "Moving #{i*b+j} (which is in fact #{i*b+jstar}) to #{i+j*a}"
    end
  end

  puts res

  # Afterwards:
  #
  #    a = 8
  #
  # 0 4  8 12 16 20 24 28
  # 2 6 10 14 18 22 26 30
  # 1 5  9 13 17 21 25 29  b = 4
  # 3 7 11 15 19 23 27 31

  # Do b a-point ffts on the rows (in-place)
  puts "Doing #{b}x #{a}-point FFT in the rows"
  b.times do |i|
    res[i*a...(i+1)*a] = fft(res[i*a...(i+1)*a])
  end

  # The sorting is afterwards:
  # because FFT maps indexes as follows:
  # 0 1 2 3 4 5 6 7 -> 0 4 2 6 1 5 3 7
  #
  # 0 16  8 24 4 20 12 28
  # 2 18 10 26 6 22 14 30
  # 1 17  9 25 5 21 13 29
  # 3 19 11 27 7 23 15 31
  #
  # Ordering of normal in-place FFT would be:
  #  0 16  8 24 4 20 12 28
  #  2 18 10 26 6 22 14 30
  #  1 17  9 25 5 21 13 29
  #  3 19 11 27 7 23 15 31
  #
  #  (which is exactly the same)

  # Voila, drink while hot.
  res
end

def ibaileys(values)
  b = values.size
  a = 1
  while a < b
    a*=2
    b/=2
  end

  # Do b a-point iffts on the rows (in-place)
  b.times do |i|
    values[i*a...(i+1)*a] = ifft(values[i*a...(i+1)*a])
  end

  res = Array.new(a*b, 0)
  omega = -2*Math::PI/values.size

  # ---
  # Multiply and transpose
  b.times do |i|
    values[i*a...(i+1)*a].each_with_index do |val, j|

      # b = 4 (index i)
      # a = 8 (index j)

      # Multiply with twiddles
      istar = bitreverse(i, b)
      k = istar*i # ??? FIXME
      angle = omega*k
      real = Math.cos(angle)
      imag = Math.sin(angle)
      t = Complex(real, imag)

      puts "Index: #{jtar},#{i}: #{val} *===* #{t} (TARGET: #{i+j*a})"

      # TODO Later: movntq
      # FIXME: Do not use jstar. Does not matter, just reorders and makes shit
      # slow (probably). Use j.
      res[i+jstar*a] = val*t
    end
  end

  # Do a b-point ffts on the rows (in-place)
  a.times do |i|
    values[i*b...(i+1)*b] = fft(values[i*b...(i+1)*b])
  end


  byebug


  # Voila, drink while the cache is still hot
  res
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

  values[0...n/2-1] = fft(values[0...n/2-1])
  values[n/2...-1]  = fft(values[n/2...-1])

  values
end

def stridedfft(values, stride)
  n = values.length/stride

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

    a = values[k*stride]
    b = values[(k+n/2)*stride]

    values[k*stride] = (a+b).to_c
    values[(k+n/2)*stride] = (a-b)*t
  end

  values[0...(n/2-1)*stride] = stridedfft(values[0...(n/2-1)*stride], stride)
  values[(n/2)*stride...-1]  = stridedfft(values[(n/2)*stride...-1], stride)

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
  a = [1, 5, 3, 6]
  stride = 4
  b = a.dup.inject([]) { |x,y| x + [y] + [0]*(stride-1)}

  puts a.to_s
  puts b.to_s

  puts "FFT:"
  a = fft(a)
  puts a
  puts "Strided:"
  b = stridedfft(b, stride)
  puts b
end

a = [1, 5, 3, 6]
b = a.dup

puts "FFT:"
a = fft(a)
puts a
puts "Baileys:"
b = baileys(b)
puts b

exit

if a*b == c
  puts "GREAT SUCCESS"
else
  puts "FAIL #{a}*#{b} = #{c} (expected: #{a*b})"
end
