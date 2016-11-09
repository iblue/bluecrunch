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

  # WARNING: Becomes unstable if alpha fluctuates!

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


puts "== ALGO 3 =="
# Subquadratic algo
def algo3(a)
  b = 10
  kt = 4 # Of k <= 4, return to basecase. kt MUST be >= 3

  # See above
  def convert_trunc(y0, k, n)
    b = 10

    # Choose FP value alpha <= log_2(b)
    alpha = Math.log(b)/Math.log(2)

    y = []
    y[0] = y0

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
      puts "  t[#{i}] = #{t[i]}"

      s[k-i] = t[i] / 2**nn(i-1, n, alpha)
      puts "  s[#{k-i}] = #{s[k-i]}"

      z[i] = t[i] % 2**nn(i-1, n, alpha)
      puts "  z[#{i}] = #{z[i]}"

      y[i] = bdiv(z[i], nn(i-1, n, alpha) - nn(i, n, alpha))
      puts "  y[#{i}] = #{y[i]}"
    end

    return s
  end

  def convert_rec(k, y, n, g)
    kt = 4 # Of k <= 4, return to basecase. kt MUST be >= 3
    b = 10

    if(k <= kt)
      return convert_trunc(y, k, n)
    else
      kh = (k+1)/2
      kl = k - kh + 1
      # Choose nh such that 4*g*b**(kh) < 2**nh
      nh = 0
      while  4*g*b**kh >= 2**nh
        nh+=1
      end

      # Choose nl such that 4*g*b**kl < 2**nl
      nl = 1
      while 4*g*b**kl >= 2**nl
        nl+=1
      end

      # yh <- floor y*2^(nh-n)
      yh = y/2**(n-nh)
      yl = bdiv(b**(k-kl)*y % 2**n, n-nl)

      sh = convert_rec(kh, yh, nh, g)
      sl = convert_rec(kl, yl, nl, g)


      def carry(s)
        # Yeah, this will be efficient in the real code
        (s.reverse.join("").to_i + 1).to_s.chars.reverse.map(&:to_i)
      end

      # if the trailing digit of sh is b-1 and the leading digit of sl is 0
      if(sh[0] == b-1 && sl[-1] == 0)
        byebug
        sh = carry(sh)
      end

      # if the trailing digit of sh is 0 and the leading digit of sl is b-1
      if(sh[0] == 0 && sl[-1] == b-1)
        byebug
        # sl <-- [0, 0, 0, ..., 0] (kl times)
        sl = [0]*kl
      end

      #??? floor(sh/b)*b^(kl) + sl
      return sl + sh[1..-1]
    end
  end

  # k = ceil(log_b(a))
  k = (Math.log(a) / Math.log(10)).ceil # Anzahl Ziffern in Basis 10
  g = [(Math.log(k)/Math.log(2)).ceil + 1, kt].max
  # Choose n such that 4*g*b^k < 2^n
  n = 0
  while 4*g*b**k >= 2**n
    n+=1
  end
  y = ((a+1)*2**n)/b**k - 1
  return convert_rec(k, y, n, g)
end

puts algo3(90090909099990909090990909090909000000999909009090900).reverse.join("")
