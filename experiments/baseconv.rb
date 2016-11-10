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
      #puts "  t[#{i}] = #{t[i]}"

      s[k-i] = t[i] / 2**nn(i-1, n, alpha)
      #puts "  s[#{k-i}] = #{s[k-i]}"

      z[i] = t[i] % 2**nn(i-1, n, alpha)
      #puts "  z[#{i}] = #{z[i]}"

      y[i] = bdiv(z[i], nn(i-1, n, alpha) - nn(i, n, alpha))
      #puts "  y[#{i}] = #{y[i]}"
    end

    return s
  end

  def convert_rec(a, k, y, n, g)
    kt = 4 # Of k <= 4, return to basecase. kt MUST be >= 3
    b = 10

    if(k <= kt)
      s = convert_trunc(y, k, n)
      return s
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

      # Not required?
      ah = (a*b**(kh - k)).floor
      al = a % b**kl

      # yh <- floor y*2^(nh-n)
      yh = y/2**(n-nh)
      yl = bdiv(b**(k-kl)*y % 2**n, n-nl)

      sh = convert_rec(ah, kh, yh, nh, g)
      sl = convert_rec(al, kl, yl, nl, g)

      def carry(s)
        add(s, [1])
      end

      def shiftdown(s)
        s[1..-1]
      end

      def shiftup(s, p)
        [0]*p + s
      end

      def add(s1,s2)
        l = [s1.length, s2.length].max

        val = (s1.reverse.join("").to_i + s2.reverse.join("").to_i)
        str = "%0#{l}d" % val
        str.chars.reverse.map(&:to_i)
      end

      def a_to_s(a, k)
        ("%0#{k}d" % a).chars.map(&:to_i).reverse
      end



      # if the trailing digit of sh is b-1 and the leading digit of sl is 0
      if(sh[0] == b-1 && sl[-1] == 0)
        sh = carry(sh)
      end

      # FIXME: Does this ever happen?
      # if the trailing digit of sh is 0 and the leading digit of sl is b-1
      if(sh[0] == 0 && sl[-1] == b-1)
        sl = [0] * kl
      end

      # Breaks if bug
      #byebug if sh != a_to_s(ah, kh)
      #byebug if sl != a_to_s(al, kl)
      # fill upper
      # sh = sh[1..-1]
      #
      ret = add(shiftup(shiftdown(sh), kl), sl)
      #byebug if ret != a_to_s(a, k)

      #??? floor(sh/b)*b^(kl) + sl
      return ret
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
  return convert_rec(a, k, y, n, g)
end

=begin
while x = rand(2**64) do
  y = algo3(x).reverse.join("").to_i
  if(x != y)
    puts x
    puts y
    puts "---"
  end
end
=end

=begin
(9999999...99999999).each do |x|
  y = algo3(x).reverse.join("").to_i
  if(x != y)
    puts x
    puts y
    puts "---"
  end
end
=end

# Subquadratic algo for fractional part
puts "== ALGO 4 =="
# Input:
# Floating point number a * 2**e (with e < 0)
# Output: Decimals
def algo4(a, e)
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
      #puts "  t[#{i}] = #{t[i]}"

      s[k-i] = t[i] / 2**nn(i-1, n, alpha)
      #puts "  s[#{k-i}] = #{s[k-i]}"

      z[i] = t[i] % 2**nn(i-1, n, alpha)
      #puts "  z[#{i}] = #{z[i]}"

      y[i] = bdiv(z[i], nn(i-1, n, alpha) - nn(i, n, alpha))
      #puts "  y[#{i}] = #{y[i]}"
    end

    return s
  end

  def convert_rec(k, y, n, g)
    kt = 4 # Of k <= 4, return to basecase. kt MUST be >= 3
    b = 10

    if(k <= kt)
      s = convert_trunc(y, k, n)
      return s
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
        add(s, [1])
      end

      def shiftdown(s)
        s[1..-1]
      end

      def shiftup(s, p)
        [0]*p + s
      end

      def add(s1,s2)
        l = [s1.length, s2.length].max

        val = (s1.reverse.join("").to_i + s2.reverse.join("").to_i)
        str = "%0#{l}d" % val
        str.chars.reverse.map(&:to_i)
      end

      def a_to_s(a, k)
        ("%0#{k}d" % a).chars.map(&:to_i).reverse
      end



      # if the trailing digit of sh is b-1 and the leading digit of sl is 0
      if(sh[0] == b-1 && sl[-1] == 0)
        sh = carry(sh)
      end

      # FIXME: Does this ever happen?
      # if the trailing digit of sh is 0 and the leading digit of sl is b-1
      if(sh[0] == 0 && sl[-1] == b-1)
        sl = [0] * kl
      end

      # Breaks if bug
      #byebug if sh != a_to_s(ah, kh)
      #byebug if sl != a_to_s(al, kl)
      # fill upper
      # sh = sh[1..-1]
      #
      ret = add(shiftup(shiftdown(sh), kl), sl)
      #byebug if ret != a_to_s(a, k)

      #??? floor(sh/b)*b^(kl) + sl
      return ret
    end
  end

  n = -e
  k = (n*Math.log(2)/Math.log(10)).ceil - 2 # Anz garantiert korrekte Ziffern in Basis 10
  g = [(Math.log(k)/Math.log(2)).ceil + 1, kt].max
  y = a

  s = convert_rec(k, y, n, g)
  puts "#{k} digits"
  return s
end

#puts algo3(2718281828459045235360287471352662497756)
#puts algo4(0x2b7e151628aed2a6abf7158809cf4f3c7, -128).reverse.join("")

puts "== ALGO 5 =="
def algo5(a, e)
  # See above
  def convert_trunc(y0, k, n)
    b = 10000_0000

    # Choose FP value alpha <= log_2(b)
    alpha = Math.log(b)/Math.log(2**32)

    y = []
    y[0] = y0

    # Write nn(i) for n - floor(i*a)
    def nn(i, n, alpha)
      n - (i*alpha).to_i
    end

    # Truncate p bits
    def bdiv(z, p)
      z/((2**32)**p)
    end

    s = []
    t = []
    z = []
    (1..k).each do |i|
      #puts "i = #{i}"
      t[i] = b*y[i-1]

      ni_1 = nn(i-1, n, alpha)
      ni = nn(i, n, alpha)

      #puts "ni_1, ni = #{ni_1}, #{ni}"
      #puts "t = #{t[i].to_s(16)}"
      #puts "  t[#{i}] = #{t[i]}"

      s[k-i] = t[i] / (2**32)**ni_1
      #puts "  s[#{k-i}] = #{s[k-i]}"

      z[i] = t[i] % (2**32)**ni_1
      #puts "  z[#{i}] = #{z[i]}"

      y[i] = bdiv(z[i], ni_1 - ni)
      #puts "  y[#{i}] = #{y[i]}"
    end

    return s
  end

  def convert_rec(k, y, n, g)
    puts "(k,y,n,g) = #{k}, #{y}, #{n}, #{g}"

    kt = 4 # Of k <= 4, return to basecase. kt MUST be >= 3
    b = 10000_0000

    if(k <= kt)
      s = convert_trunc(y, k, n)
      return s
    else
      kh = (k+1)/2
      kl = k - kh + 1
      # Choose nh such that 4*g*b**(kh) < 2**nh
      nh = 0
      while  4*g*b**kh >= (2**32)**nh
        nh+=1
      end

      # Choose nl such that 4*g*b**kl < 2**nl
      nl = 1
      while 4*g*b**kl >= (2**32)**nl
        nl+=1
      end

      # Truncate p bits
      def bdiv(z, p)
        z/((2**32)**p)
      end

      # yh <- floor y*2^(nh-n)
      yh = (y/(2**32)**(n-nh)).to_i
      yl = (bdiv(b**(k-kl)*y % (2**32)**n, n-nl)).to_i

      byebug
      sh = convert_rec(kh, yh, nh, g)
      sl = convert_rec(kl, yl, nl, g)

      def carry(s)
        add(s, [1])
      end

      def shiftdown(s)
        s[1..-1]
      end

      def shiftup(s, p)
        [0]*p + s
      end

      def add(s1,s2)
        b = 10000_0000
        l = [s1.length, s2.length].max

        carry = 0
        s = []
        (0..l).each do |i|
          sum = (s1[i] || 0) + (s2[i] || 0) + carry
          carry = sum / b
          s[i]  = sum % b
        end

        s
      end

      def a_to_s(a, k)
        ("%0#{k}d" % a).chars.map(&:to_i).reverse
      end

      # if the trailing digit of sh is b-1 and the leading digit of sl is 0
      if(sh[0] == b-1 && sl[-1] == 0)
        sh = carry(sh)
      end

      # FIXME: Does this ever happen?
      # if the trailing digit of sh is 0 and the leading digit of sl is b-1
      if(sh[0] == 0 && sl[-1] == b-1)
        sl = [0] * kl
      end

      # Breaks if bug
      #byebug if sh != a_to_s(ah, kh)
      #byebug if sl != a_to_s(al, kl)
      # fill upper
      # sh = sh[1..-1]
      #
      ret = add(shiftup(shiftdown(sh), kl), sl)
      #byebug if ret != a_to_s(a, k)

      #??? floor(sh/b)*b^(kl) + sl
      return ret
    end
  end

  #n = -e
  #k = (n*Math.log(2**32)/Math.log(100000000)).ceil - 2 # Anz garantiert korrekte Ziffern in Basis 10
  #y = a

  n = 22
  k = 25
  g = 4
  y = a

  s = convert_rec(k, y, n, g)
  puts "#{k} digits"
  return s
end

s = algo5(228776509382463697450290884140869401551305345104850785467877452727290000022059215493424723739701291948663852630304610978205778895706868485166037611507017266201971612553383158879956605321441860961188689115602691286, -22)
puts s.map{|d| "%08d" % d}.reverse.join("")
