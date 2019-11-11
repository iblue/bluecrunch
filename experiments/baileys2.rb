# Baileys Algorithmn Prototype

require "byebug"

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

  # r 25441541210108625536
  # e 2545553108625536

  # 27
  #  45
  #   54
  #    27
  #     9

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
  a = fft(a)
  b = fft(b)
  c = fft_pointwise(a,b)
  c = ifft(c)
  c = fft_to_int(c)
end

a = 1239912
b = 9992332
c = mul(a,b)

if a*b == c
  puts "GREAT SUCCESS"
else
  puts "FAIL #{a}*#{b} = #{c} (expected: #{a*b})"
end
