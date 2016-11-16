# Baileys Algorithmn Prototype

require "symbolic"
require "byebug"

values = (0..15).map{|n| "i#{n}"}.map{|name| var :name => name}

def fft(values)
  n = values.length

  if n == 1
    return values
  end

  omega = var :name => "ω"

  (0..n/2-1).each do |k|
    a = values[k]
    b = values[k+n/2]

    values[k] = a+b
    values[k+n/2] = (a-b)*(omega**(k.to_r/n))
  end

  values[0..n/2-1] = fft(values[0..n/2-1])
  values[n/2..-1]  = fft(values[n/2..-1])

  values
end

#puts fft(values)

def baileys(values)
  #[
  #  [a,b,c,d],
  #  [e,f,g,h]
  #]

  n1 = values[0].size
  n2 = values.size

  # Do FFT on the cols
  (0..n1-1).each do |col_idx|
    # Do multiple in parallel to better use cache
    col = values.map{|x| x[col_idx]}
    col = fft(col)

    puts col.inspect

    (0..n2-1).each do |row_idx|
      values[row_idx][col_idx] = col[row_idx]
    end
  end

  # Multiply by twiddles
  omega = var :name => "ω"
  (0..n1-1).each do |col_idx|
    (0..n2-1).each do |row_idx|
      values[row_idx][col_idx] *= omega**(col_idx*row_idx.to_r/(n1*n2))
    end
  end

  # Transpose
  #values = values.transpose
  # => Fuck this. We are doing fft on the rows which is faster anyway

  values = values.map { |row| fft(row) }

  values
end

matrix = [
  [var(:name => 'i0'), var(:name => 'i1'), var(:name => 'i2'), var(:name => 'i3')],
  [var(:name => 'i4'), var(:name => 'i5'), var(:name => 'i6'), var(:name => 'i7')],
  [var(:name => 'i8'), var(:name => 'i9'), var(:name => 'i10'), var(:name => 'i11')],
  [var(:name => 'i12'), var(:name => 'i13'), var(:name => 'i14'), var(:name => 'i15')],
]

puts "== Fenske-Baileys =="
puts baileys(matrix)

puts "== Normal FFT =="
puts fft(values)
