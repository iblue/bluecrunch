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


a = [4, 3, 2, 1, 0, 0, 0, 0]
b = [3, 2, 1, 0, 0, 0, 0, 0]

puts a
puts b
puts "---"
a = fft(a)
puts a
b = fft(b)
puts b
puts "---"
c = a.zip(b).map{|x,y| x*y}
c = ifft(c)
c = c.map{|x| x*0.125}

# 1
#  4
#  10
#   16
#    17
#     12
#= 151782

puts c
