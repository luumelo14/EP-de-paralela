#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define X 50
#define Y 100
#define Z 12.4
#define NITER 1000

float** criaMatriz (int m, int n);
int temGota (int maxRand);
int inteiroAleatorio (int maximo);

typedef struct {
	int x;
	int y;
	int tempo;
} gota;


int main (int argc, char* argv[]) {

	int i, numGotas = 0;
	gota * gotas;

	gotas = malloc (NITER * sizeof(gota));
	criaMatriz (X, Y);
	
	for(i = 0; i < NITER; i++) {
		if(rand() <= (Z/100) * RAND_MAX) {
			gotas[numGotas].x = inteiroAleatorio(X);

			printf("%d \n", gotas[numGotas].x);
		}
		else {
			//segundo de sol
		}
	}
	
	return 0;
}

float** criaMatriz (int m, int n)
{
    float** M; 
	int i, j;
	M = malloc (m * sizeof (float *));
	for (i = 0; i < m; ++i)
	   M[i] = calloc (n, sizeof (float));

	return M;
}


int inteiroAleatorio (int maximo)
{
    return ( (double) rand() / (double) (RAND_MAX + 1)) * maximo;
}