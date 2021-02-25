# Source directory
SRC_DIR := ../../src/Serial

INC := -I$(SRC_DIR)
LIB := -lm

CC	:= gcc-10
# CC  := clang
# CC	:= icc

CFLAGS := -O3 -Winline -Wall
# CFLAGS := $(CFLAGS) -march=native

# Flag to turn off the correction step via CDF evaluation
# CFLAGS := $(CFLAGS) -D NO_CDF

# Flag to specify whether to run binoinv_timing with just DOUBLE or SINGLE routines
# CFLAGS := $(CFLAGS) -D DOUBLE
# CFLAGS := $(CFLAGS) -D SINGLE

# Flag to turn on the checks on the binocinv function
CFLAGS := $(CFLAGS) -D COMPLEMENTARY

# Flag to turn on the checks on the vector versions
CFLAGS := $(CFLAGS) -D VECTOR

# Flag to turn off the checks on the scalar versions
# CFLAGS := $(CFLAGS) -D SKIP_SCALAR

# Flag to speed up the binoinv_check code by using shortcuts
CFLAGS := $(CFLAGS) -D FAST_CHECK
ifeq ($(CC), icc)
	# When using FAST_CHECK test need to tell Intel to not flush to zero
	CFLAGS := $(CFLAGS) -no-ftz
endif

# Flags to test accuracy by increments on a linear scale, rather than the standard
# increments on a log scale
# CFLAGS := $(CFLAGS) -D PLUS
# CFLAGS := $(CFLAGS) -D COUNT_P=100000

# Dependencies
DEPS := $(SRC_DIR)/binoinv.h \
		$(SRC_DIR)/binocdf.h \
		$(SRC_DIR)/binoinvf.h \
		$(SRC_DIR)/binocdff.h

all:	binoinv_timing binoinv_check

binoinv_timing:	    binoinv_timing.c $(DEPS) Makefile
	${CC} $(CFLAGS) binoinv_timing.c -o binoinv_timing $(INC) $(LIB)

binoinv_check:		binoinv_check.c binoinv_quad.o $(DEPS) Makefile
	${CC} $(CFLAGS) binoinv_check.c binoinv_quad.o -o binoinv_check $(INC) $(LIB) -lquadmath

binoinv_quad.o:	binoinv_quad.c Makefile
# CLANG does not find libquadmath so compile the poissinv_quad with gcc instead
ifeq ($(CC), clang)
	gcc $(CFLAGS) -c binoinv_quad.c
else
	${CC} $(CFLAGS) -c binoinv_quad.c
endif

# Timing routine to check performance across different compilers
timing_compilers:   binoinv_timing.c $(DEPS) Makefile
	gcc-10 $(CFLAGS) binoinv_timing.c -o binoinv_timing_gcc   $(INC) $(LIB)
	clang  $(CFLAGS) binoinv_timing.c -o binoinv_timing_clang $(INC) $(LIB)
	icc    $(CFLAGS) binoinv_timing.c -o binoinv_timing_icc   $(INC) $(LIB)

# Test routine to check whether exponentiation of log or repeated mulitplication is faster
pow_test: 			pow_test.c Makefile
	${CC} $(CFLAGS) pow_test.c -o pow_test $(LIB)

clean:
	rm -f binoinv_timing binoinv_check binoinv_quad.o pow_test binoinv_timing_*
