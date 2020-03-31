CC=gcc
CFLAGS=-Wall -O3 -mfma -frounding-math -ffloat-store -march=native
CLIBS=-lgmp -lmpfr -lm
CMAKE_PATH=
CMAKE_INCLUDE=
DEPS=stochrndhp.h stochrnddw.h
OBJ=stochrndhp.o stochrnddw.o

all: test_performance test_validation

run_%: test_%
	./$<

test_%: test_%.o $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(CLIBS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm *.o test_performance test_validation
