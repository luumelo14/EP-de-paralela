CC = gcc
CFLAGS = -Wall -fopenmp
LDLIBS = -lm -fopenmp

ondas: ondas.o

ondasantigo: ondasantigo.o

ondasintermediario: ondasintermediario.o

.PHONY clean:
	rm *.o ondas ondasantigo ondasintermediario
