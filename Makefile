MPICXX = mpicxx
RM = rm -f
PROJECT_ROOT := .

CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -std=c++11 -fopenmp
CFLAGS =	-O2 -g -Wall -fmessage-length=0 -std=c99 -funroll-loops
LDFLAGS = -lpthread -lrt -lfftw3 -lgmp -lgmpxx -lm 
MPIFLAGS = -I$(PROJECT_ROOT)/include -L$(PROJECT_ROOT)/lib \
  -lboost_mpi -lboost_serialization -lboost_program_options
SRCS = mocksig.cpp fft.cpp reduction.cpp qdsawrapper.cpp
OBJS = $(subst .cpp,.o,$(SRCS))
FE25519OBJ = obj/consts.o \
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

C25519 = $(FE25519OBJ) \
		 obj/ladder.o \
		 obj/dh.o \
		 obj/sign.o \
		 obj/hash.o

FFT_BIN = test_fft
ATTACK_MPI_BIN = attack_mpi
SIGGEN_MPI_BIN = siggen_mpi

ALL_TARGETS = $(FFT_BIN) $(ATTACK_MPI_BIN) $(SIGGEN_MPI_BIN)

$(FFT_BIN): mocksig.o fft.o test_fft.cpp obj/c25519.a
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	
$(ATTACK_MPI_BIN): mocksig.o fft.o reduction.o attack_mpi.cpp obj/c25519.a
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPIFLAGS)

$(SIGGEN_MPI_BIN): mocksig.o qdsawrapper.o siggen_mpi.cpp obj/c25519.a
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(MPIFLAGS)

mocksig.o: mocksig.h
fft.o: mocksig.h fft.h
reduction.o: reduction.cpp
	$(MPICXX) $(CXXFLAGS) -c -o $@ $^ $(LDFLAGS) $(MPIFLAGS)
qdsawrapper.o: qdsawrapper.h obj/c25519.a
obj/%.o: qDSA/Curve25519-asm/%.[csS]
	mkdir -p obj/
	$(CC) -fPIC $(CFLAGS) -c -o $@ $^
obj/c25519.a: $(C25519)
	$(AR) -ar cr $@ $^


.PHONY: all clean

all: $(ALL_TARGETS)

clean:
	$(RM) $(OBJS) $(ALL_TARGETS) $(C25519) $(FE25519OBJ)
