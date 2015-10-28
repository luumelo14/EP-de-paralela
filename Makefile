CC = gcc
CFLAGS = -Wall -ansi -pedantic
LDLIBS = -lm

ondas: ondas.o

.PHONY clean:
	rm *.o ondas
