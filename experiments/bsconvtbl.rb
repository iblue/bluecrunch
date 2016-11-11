def conv(k)
  return if k < 3

  kh = (k+1)/2;
  kl = k - kh + 1;

  $a << k-kl
  conv(kl)
  conv(kh) if kl != kh
end

def bconv(k)
  [23]
end

(3..100).each do |k|
  $a = []
  conv(k)
  aa = $a.uniq.sort
  puts "#{k.to_s(2)} -> #{aa.inspect} [#{((k-1)).to_s(2)}]"
  bb = bconv(k).sort
  puts "#{k.to_s(2)} -> #{bb.inspect}"
end
