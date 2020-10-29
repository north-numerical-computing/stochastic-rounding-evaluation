CC=gcc-mp-10
# We now use "-mfpmath=sse -msse2" to avoid x87 80-bit arithmetic.
WCFLAGS=-Wall -Wextra # Warnings
MCFLAGS=-march=native -mfpmath=sse # -mfma -msse2 included in -march=native.
CFLAGS=$(WCFLAGS) $(MCFLAGS)
CFLAGS_LO=$(CFLAGS) -O0 # Use with switching of rounding modes.
CFLAGS_HO=$(CFLAGS) -O3 # Use with RN-only and MPFR version of algorithms.
CLIBS=-lgmp -lmpfr -lm
CMAKE_PATH=/opt/local/lib/
CMAKE_INCLUDE=/opt/local/include/
DEPS=stochrndhp.h stochrnddw.h
OBJ=stochrndhp.o stochrnddw_rn.o stochrnddw_ar.o
.PRECIOUS: %.o




all: test_performance test_validation

run_%: test_%
	./$<

test_%: test_%.o $(OBJ)
	$(CC) $(CFLAGS_HO) -L$(CMAKE_PATH) -o $@ $^ $(CLIBS)

stochrnddw_ar.o: stochrnddw_ar.c stochrnddw.h
	$(CC) $(CFLAGS_LO) -I$(CMAKE_INCLUDE) -c -o $@ $<

stochrnddw_rn.o: stochrnddw_rn.c stochrnddw.h
	$(CC) $(CFLAGS_HO) -I$(CMAKE_INCLUDE) -c -o $@ $<

stochrndhp.o: stochrndhp.c stochrndhp.h
	$(CC) $(CFLAGS_HO) -I$(CMAKE_INCLUDE) -c -o $@ $<

test_%.o: test_%.c $(DEPS)
	$(CC) $(CFLAGS_HO) -I$(CMAKE_INCLUDE) -c -o $@ $<




profile: test_performance
	valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes \
		./test_performance
	qcachegrind callgrind.out.*

query:
	gcc -Q -march=native --help=target

clean:
	rm -f *.o test_performance test_validation callgrind.out.*
