require "byebug"

# Load {bits} bits from {word} word
def get(bits, word)
  puts "get #{bits} from #{word}"
end

# Load word if exists
def ld(word)
  puts "load #{word}"
end

def conv(pwcnt, pwpos, cwcnt, cwpos)
  if pwcnt == cwcnt
    get(cwpos-pwpos, cwcnt)
  else
    bits = 32*cwcnt + cwpos - 32*pwcnt - pwpos

    get(32-pwpos, pwcnt)
    ld(cwcnt) if cwpos%bits != 0
    get(cwpos%bits, cwcnt) if cwpos%bits != 0
  end
end

def gencode(bits_per_point)
  # get 12 bit from word 0 => 0 (20 left)
  # get 12 bit from word 0 => 1 ( 8 left)
  # get  8 bit from word 0 => 2 ( 0 left)
  # get  4 bit from word 1 => 2 (
  # get 12 bit from word 1 => 3
  # get 12 bit from word 1

  bitpos  = 0
  while true do
    pwcnt    = bitpos/32
    pwpos    = bitpos%32
    bitpos += bits_per_point
    cwcnt    = bitpos/32
    cwpos    = bitpos%32

    #puts "#{bitpos}: (#{pwcnt},#{pwpos}) (#{cwcnt},#{cwpos})"
    conv(pwcnt, pwpos, cwcnt, cwpos)

    if cwpos == 0
      break
    end
  end
end

gencode(9)
