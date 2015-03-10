#CC=icc
CC=gcc
#CC=clang
VALGRIND=valgrind --tool=memcheck --leak-check=yes
STRIP=strip
CFLAGS=-std=c11 -Wall -Werror -ggdb -O3 -fopenmp -msse3
BINARY=bluecrunch
LIBS=-lm
ifeq ($(CC),icc)
	CFLAGS += -ipo
endif

INCLUDES = -I./include -I.
SOURCES  = $(shell find -path "./bigfloat/*" -name "*.c")
SOURCES += $(shell find -path "./fft/*"      -name "*.c")
SOURCES += main.c
HEADERS  = $(shell find -path "./include/*"  -name "*.h")
OBJECTS  = $(SOURCES:.c=.o)

$(BINARY): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(BINARY) $(OBJECTS) $(LIBS)

.c.o: $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(BINARY) $(OBJECTS)

.PHONY: stats
stats:
	git ls-files | xargs wc -l

.PHONY: build
build: CFLAGS += -DNDEBUG
build: clean $(BINARY)

.PHONY: debug
debug: CFLAGS := $(filter-out -O3,$(CFLAGS)) -DDEBUG -O0
ifeq ($(CC),gcc)
debug: CFLAGS += -fsanitize=address
endif
debug: clean $(BINARY)

# Build static, pack and remove traces of UPX
release: CFLAGS=-static -std=c11 -Wall -Werror -O3 -fopenmp -msse3
# static doesn't work with icc?
release: CC=gcc
release: clean build
	$(STRIP) --strip-all $(BINARY)
	strip --remove-section=.note.gnu.build-id $(BINARY)
	upx --best $(BINARY)
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

