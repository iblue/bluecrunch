# Baileys Algorithmn Prototype

require "symbolic"
require "byebug"

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

# Bit reverse ordering
def reorder(arr)
  raise if arr.size != 4
  out = [arr[0],arr[2],arr[1],arr[3]]
end

def baileys(values)
  #[
  #  [a,b,c,d],
  #  [e,f,g,h]
  #]

  n1 = values[0].size
  n2 = values.size

  # Do FFT on the cols
  values = values.transpose
  values = values.map { |row| reorder(fft(row)) }
  values = values.transpose

  # Multiply by twiddles
  omega = var :name => "ω"
  (0..n1-1).each do |col_idx|
    (0..n2-1).each do |row_idx|
      values[row_idx][col_idx] *= omega**(col_idx*row_idx.to_r/(n1*n2))
    end
  end

  # Do FFT on the rows
  values = values.map { |row| fft(row) }

  # Skip re-transpose, because we do not care about ordering anyway.
  values
end

if true
  matrix = [
    [var(:name => 'i0'), var(:name => 'i1'), var(:name => 'i2'), var(:name => 'i3')],
    [var(:name => 'i4'), var(:name => 'i5'), var(:name => 'i6'), var(:name => 'i7')],
    [var(:name => 'i8'), var(:name => 'i9'), var(:name => 'i10'), var(:name => 'i11')],
    [var(:name => 'i12'), var(:name => 'i13'), var(:name => 'i14'), var(:name => 'i15')],
  ]

  values = (0..15).map{|n| "i#{n}"}.map{|name| var :name => name}
else
  matrix = [
    [var(:name => 'i0'), var(:name => 'i1')],
    [var(:name => 'i2'), var(:name => 'i3')],
    [var(:name => 'i4'), var(:name => 'i5')],
    [var(:name => 'i6'), var(:name => 'i7')],
  ]

  values = (0..7).map{|n| "i#{n}"}.map{|name| var :name => name}
end

puts "== Fenske-Baileys =="
puts baileys(matrix)

puts "== Normal FFT =="
puts fft(values)
