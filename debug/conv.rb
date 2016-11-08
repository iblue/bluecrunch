def to_base(value, base, rep)
  i = 0;
  while(value > 0) do
    if rep == 16
      puts "a->coef[%d] = 0x%08x;" % [i, value % base]
    else
      puts "a->coef[%d] = %d;" % [i, value % base]
    end
    value /= base
    i+=1
  end
  nil
end

def to16(value)
  to_base(value, 2**32, 16)
end

def to10(value)
  to_base(value, 100000000, 10)
end
