#CC=/opt/intel/bin/icc
#CFLAGS=-std=c11 -Wall -Werror -ggdb -O3 -mavx
#CC=gcc-5 # Will run into inling errors with cilk plus and intrinsics
#CC=/opt/magic/bin/clang # Needs to be compiled with cilk support (see https://cilkplus.github.io/)
#CC=/opt/gcc-trunk/bin/gcc
#CFLAGS=-std=c11 -Wall -Werror -ggdb -O3 -mavx -fcilkplus -fprofile-use
CC=gcc
CFLAGS=-std=c11 -Wall -Werror -ggdb -O3 -mavx -fcilkplus -funroll-loops -fprofile-use=. -flto
VALGRIND=valgrind --tool=memcheck --leak-check=yes
STRIP=strip
BINARY=bluecrunch
LIBS=-lm #-ltcmalloc # Needs libgoogle-perftools-dev
ifeq ($(CC),icc)
	CFLAGS += -ipo
endif
ifeq ($(CC),gcc)
	# => Infinite loop in calc of number of terms
	#CFLAGS += -ffast-math
endif

MAIN_C = main.c
MAIN_O = $(MAIN_C:.c=.o)

BENCH_C = bench.c
BENCH_O = $(BENCH_C:.c=.o)

INCLUDES = -I./include -I.
SOURCES  = $(shell find -path "./bigfloat/*" -name "*.c")
SOURCES += $(shell find -path "./fft/*"      -name "*.c")
SOURCES += $(shell find -path "./util/*"     -name "*.c")
HEADERS  = $(shell find -path "./include/*"  -name "*.h")
OBJECTS  = $(SOURCES:.c=.o)

TEST_SOURCES = $(shell find -path "./test/*"     -name "*.c")
TEST_OBJECTS = $(TEST_SOURCES:.c=.o)
TESTS        = $(TEST_SOURCES:.c=.test)

.PRECIOUS: %.c %.o %.h

$(BINARY): codegen $(OBJECTS) $(MAIN_O)
	$(CC) $(CFLAGS) -o $(BINARY) $(OBJECTS) $(MAIN_O) $(LIBS)

bench: $(OBJECTS) $(BENCH_O)
	$(CC) $(CFLAGS) -o bench $(OBJECTS) $(BENCH_O) $(LIBS)

codegen:
	make -C experiments/codegen

.c.o: $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(BINARY) bench $(OBJECTS) $(TEST_OBJECTS) $(TESTS) $(MAIN_O)
	make -C experiments/codegen clean

.PHONY: stats
stats:
	git ls-files | xargs wc -l

.PHONY: build
build: CFLAGS += -DNDEBUG
build: clean $(BINARY)

.PHONY: tests run_tests
test: $(TESTS) run_tests

run_tests: $(TESTS)
	$(foreach test,$(TESTS),echo $(test); $(test);)

%.test: %.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $< $(OBJECTS) $(INCLUDES) $(LIBS)

.PHONY: debug
debug: CFLAGS := $(filter-out -O3,$(CFLAGS)) -DDEBUG -O0
ifeq ($(CC),gcc)
debug: CFLAGS += -fsanitize=address
endif
debug: $(BINARY) bench

.PHONY: debug-test
debug-test: CFLAGS := $(filter-out -O3,$(CFLAGS)) -DDEBUG -O0
ifeq ($(CC),gcc)
debug-test: CFLAGS += -fsanitize=address
endif
debug-test: test

# Build static, pack and remove traces of UPX
release: CFLAGS=-static-intel -std=c11 -Wall -Werror -O3 -mavx
# static doesn't work with icc?
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

