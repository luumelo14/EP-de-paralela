CC = gcc
CFLAGS = -Wall -fopenmp -O3
LDLIBS = -lm -fopenmp

ondas: ondas.o

ondasantigo: ondasantigo.o

ondasintermediario: ondasintermediario.o

ondasfinal: ondasfinal.o 

.PHONY clean:
	rm *.o  ondasfinal
