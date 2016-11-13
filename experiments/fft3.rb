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

puts len(11)
puts len(33)
