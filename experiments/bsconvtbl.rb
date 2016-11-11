def conv(k)
  return if k < 3

  kh = (k+1)/2;
  kl = k - kh + 1;
  exp = k - kl

  conv(kh)
  conv(kl) if kh != kl

  unless $a.include?(exp)
    $a << exp
  end

end

(3..100).each do |k|
  $a = []
  conv(k)
  aa = $a
  puts "#{k} -> #{aa.inspect.sort}"
end
