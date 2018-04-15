CC=gcc
CFLAGS  = -g -O3 -Wall -Wextra -Werror

all: test/test

fe25519obj = obj/consts.o \
	obj/fe25519_add.o \
	obj/fe25519_copy.o \
	obj/fe25519_freeze.o \
	obj/fe25519_getparity.o \
	obj/fe25519_invert.o \
	obj/fe25519_iseq.o \
	obj/fe25519_iszero.o \
	obj/fe25519_mul.o \
	obj/fe25519_neg.o \
	obj/fe25519_pack.o \
	obj/fe25519_pow2523.o \
	obj/fe25519_setint.o \
	obj/fe25519_square.o \
	obj/fe25519_sub.o \
	obj/fe25519_unpack.o \
	obj/sc25519_add.o \
	obj/sc25519_barrett.o \
	obj/sc25519_from32bytes.o \
	obj/sc25519_from64bytes.o \
	obj/sc25519_from_shortsc.o \
	obj/sc25519_iszero.o \
	obj/sc25519_lt.o \
	obj/sc25519_mul.o \
	obj/sc25519_mul_shortsc.o \
	obj/sc25519_slide.o \
	obj/sc25519_sub_nored.o \
	obj/sc25519_to32bytes.o \
	obj/sc25519_window4.o \
	obj/ull4_mul.o \
	obj/keccak.o \
	obj/KeccakSponge.o \
	obj/KeccakF-1600-x86-64-asm.o

c25519 = $(fe25519obj) \
				 obj/ladder.o \
				 obj/dh.o \
				 obj/sign.o \
				 obj/hash.o

test/test: test/test.c test/print.c obj/c25519.a
	$(CC) $(CFLAGS) $^ -o $@

test/print: test/print.c obj/c25519.a
	$(CC) $(CFLAGS) $^ -o $@

obj/c25519.a: $(c25519)
	$(AR) -ar cr $@ $^

obj/%.o: %.[cSs]
	mkdir -p obj/
	$(CC) -fPIC $(CFLAGS) -c $^ -o $@

.PHONY: clean

clean:
	-rm -r obj/*
	-rm -r test/test

run:
	test/test
