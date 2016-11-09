# Implements stuff from
# Cyril Bouvier, Paul Zimmermann. Division-Free Binary-to-Decimal Conversion, 2013


# Algorithmn 1
def algo1(a, b, k)
  # Choose an integer n such that 2*b**k < 2**n
  n = 0
  while 2*b**k >= 2**n
    n+=1
  end
  y = []
  y[0] = ((a+1)*2**n)/b**k - 1
  puts "y[0] = #{y[0]}"

  s = []
  t = []
  (1..k).each do |i|
    t[i] = b*y[i-1]
    puts "t[#{i}] = #{t[i]}"
    s[k-i] = t[i]/2**n
    puts "s[#{k-i}] = #{s[k-i]}"
    y[i] = t[i]%2**n
    puts "y[#{i}] = #{y[i]}"
  end
  return s
end

puts algo1(271828, 10, 8).join(",")

