# Baileys Algorithmn Prototype

require "symbolic"
require "byebug"

values = ('a'..'d').map{|name| var :name => name}

def fft(values)
  n = values.length

  if n == 1
    return values
  end

  omega = var :name => "ω_#{n}"

  (0..n/2-1).each do |k|
    a = values[k]
    b = values[k+n/2]

    values[k] = a+b
    values[k+n/2] = (a-b)*(omega**k)
  end

  values[0..n/2-1] = fft(values[0..n/2-1])
  values[n/2..-1]  = fft(values[n/2..-1])

  values
end

puts fft(values)
