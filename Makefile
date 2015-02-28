CPP=icc
STRIP=strip
CFLAGS=-std=c++11 -Wall -Werror -msse3 -O3 -fopenmp -ipo -mtune=native -march=native
BINARY=bluecrunch
LIBS=-lm

$(BINARY): main.o fft.o bigfloat.o
	$(CPP) $(CFLAGS) -o $(BINARY) main.o fft.o bigfloat.o $(LIBS)

.cpp.o: $(HEADERS)
	$(CPP) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(BINARY) main.o fft.o bigfloat.o

.PHONY: stats
stats:
	git ls-files | xargs wc -l

.PHONY: build
build: CFLAGS += -DNDEBUG
build: clean $(BINARY)

# Build static, pack and remove traces of UPX
release: CFLAGS += -static
release: build
	$(STRIP) --strip-all $(BINARY)
	strip --remove-section=.note.gnu.build-id $(BINARY)
	upx --ultra-brute --best $(BINARY)
	sed -i 's/UPX!/B0Un/g' $(BINARY)
	sed -i 's/PROT_EXEC/\xff\xff\xff\x00\xff\xff\x9a\xff\x00/g' $(BINARY)
	sed -i 's/PROT_WRITE fail/\xff\xff\xff\x00\xff\xff\xff\xff\xff\x09 6@AP/g' $(BINARY)
	sed -i 's/This file is packed/\xff\x00\xff\xff\xff\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff/g' $(BINARY)
	sed -i 's/Info: /\xff\xff\xff\xff\xff\xff/g' $(BINARY)
	sed -i 's/with the UPX exec/\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff exec/g' $(BINARY)
	sed -i 's/utable packer http:/\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff/g' $(BINARY)
	sed -i 's/upx\.sf\.net/\xff\xff\xffB\xff\xffZ\xff\xff\xff/g' $(BINARY)
	sed -i 's/UPX 3\.91 Copyright/\xff\xff\xff \xff\xff\xff\xff \xff\xff\xff\xff\xff\xff\xff\xff\xff/g' $(BINARY)
	sed -i 's/1996-2013 the UPX Team/\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff/g' $(BINARY)
	sed -i 's/All Rights Reserved/4ll\x00r1tes\x00\x00r\x03vers\x03d/g' $(BINARY)
	sed -i 's/GCC: (Ubuntu 4/iblue pwnz :)\x00/g' $(BINARY)

