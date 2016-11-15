# Baileys Algorithmn Prototype

require "symbolic"

omega = var :name => 'ω'
a = var :name => 'a'
b = var :name => 'b'
c = var :name => 'c'
d = var :name => 'd'

at = a +            b +            c +            d
bt = a + omega    * b + omega**2 * c + omega**3 * d
ct = a + omega**2 * b +            c + omega**2 * d
dt = a + omega**3 * b + omega**2 * c + omega**1 * d

puts at, bt, ct, dt
