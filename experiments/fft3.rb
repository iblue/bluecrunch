def len(n)
  k = 1
  while true
    if 2*k > n
      return 2*k
    end
    if 3*k > n
      return 3*k
    end
    if 5*k > n
      return 5*k
    end
    k *= 2
  end
end

i = 0
(0..63).each do |k|
  puts "#{2*2**k}, // #{i} => 2*2^#{k}"
  i += 1
  puts "#{3*2**k}, // #{i} => 3*2^#{k}"
  i += 1
end
