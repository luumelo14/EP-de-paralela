CC = gcc
CFLAGS = -Wall -fopenmp
LDLIBS = -lm -fopenmp

ondas: ondas.o

.PHONY clean:
	rm *.o ondas
