CC = gcc
CFLAGS = -Wall
LDLIBS = -lm

ondas: ondas.o

.PHONY clean:
	rm *.o ondas
