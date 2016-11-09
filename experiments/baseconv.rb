# Implements stuff from
# Cyril Bouvier, Paul Zimmermann. Division-Free Binary-to-Decimal Conversion, 2013
require "byebug"

puts "== ALGO 1 =="

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

puts "== ALGO 2 =="

# Algorithmn 2
def algo2(a, b, k)
  # Choose an integer n such that 2*b**k < 2**n
  n = 0
  while 2*b**k >= 2**n
    n+=1
  end
  y = []
  y[0] = ((a+1)*2**n)/b**k - 1
  puts "y[0] = #{y[0]}"

  # Choose FP value alpha <= log_2(b)
  alpha = Math.log(b)/Math.log(2)

  # Write nn(i) for n - floor(i*a)
  def nn(i, n, alpha)
    n - (i*alpha).to_i
  end

  # Truncate p bits
  def bdiv(z, p)
    z/2**p
  end

  s = []
  t = []
  z = []
  (1..k).each do |i|
    t[i] = b*y[i-1]
    puts "t[#{i}] = #{t[i]}"

    s[k-i] = t[i] / 2**nn(i-1, n, alpha)
    puts "s[#{k-i}] = #{s[k-i]}"

    z[i] = t[i] % 2**nn(i-1, n, alpha)
    puts "z[#{i}] = #{z[i]}"

    y[i] = bdiv(z[i], nn(i-1, n, alpha) - nn(i, n, alpha))
    puts "y[#{i}] = #{y[i]}"
  end
  return s
end

puts algo2(271828, 10, 7).join(",")
