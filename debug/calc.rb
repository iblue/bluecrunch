File.open('out.test').each do |line|
  a,b,c = line.split(",")
  a = a.to_i
  b = b.to_i
  c = c.to_i
  if(a.to_i*b.to_i != c.to_i)
    puts "    #{a.to_i} * #{b.to_i} != #{c.to_i}"
    puts "<=> #{a} * #{b} != #{c}"
  end
end
