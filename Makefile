CC = gcc
CFLAGS = -Wall -fopenmp
LDLIBS = -lm -fopenmp

ondas: ondas.o

ondasantigo: ondasantigo.o

.PHONY clean:
	rm *.o ondas ondasantigo
