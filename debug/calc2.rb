@mode = nil

File.open('debug-20150404.log').each do |line|
  if line =~ /^mul/
    @mode = :mul
    @line_mode = line
    next
  end

  if line =~ /^add/
    @mode = :add
    @line_mode = line
    next
  end

  if line =~ /coef = { (.+) }/
    @coef = $1.split(", ")
    @num  = @coef.map do |value|
      value.to_i(16)
    end.each_with_index.map do |value, index|
      value*(2**(32*index))
    end.reduce(&:+)
  end

  if line =~ /^a =/
    @a = @num
    @line_a = line
    next
  end

  if line =~ /^b =/
    @b = @num
    @line_b = line
    next
  end

  if line =~ /^t =/
    @t = @num
    @line_t = line
  end

  case @mode
  when :mul
    @expected = @a*@b
  when :add
    @expected = @a+@b
  else
    raise
  end

  if(@t != @expected)
    puts "Got 0x%x but expected 0x%x" % [@expected, @t]
    puts @line_mode
    puts @line_a
    puts @line_b
    puts @line_t
    break
  end
end
