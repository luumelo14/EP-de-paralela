#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* #define H 5
// #define L 5
// #define P 90
// #define NITER 30
// #define T 120
// #define V 10
// #define ALT 5
// #define LARG 5 */


float** criaMatriz (int m, int n);
int temGota (int maxRand);
int inteiroAleatorio (int maximo);
float calculaDistancia (int x1, int y1, int x2, int y2, float aspectx, float aspecty);
void escreveArquivo(float** lago);
void pegaEntrada(int argc, char* argv[]);

typedef struct {
	int x;
	int y;
	float tempo;
} gota;

int H, L, NITER, T, SEED, NPROCS;
float P, V, ALT, LARG, EPS;

int main (int argc, char* argv[]) {

	int i, j, k, l, numGotas = 0;
	float t, aux;
	gota * gotas;
	float** lago;
	float aspectx, aspecty;

	pegaEntrada(argc, argv);
	aspectx = (float) LARG/L;
	aspecty = (float) ALT/H;
	gotas = malloc (NITER * sizeof(gota));
	lago = criaMatriz (H, L);
	
	for(i = 0; i < NITER; i++) {
		if(numGotas != 0) {
			for(j = 0; j < H; j++) {
				for(k = 0; k < L; k++) {
					lago[j][k] = 0;
					for(l = 0; l < numGotas; l++) {
						t = (float) i * T / NITER - gotas[l].tempo;
						/* aux = ro - v*t */
						aux = calculaDistancia(j, k, gotas[l].x, gotas[l].y, aspectx, aspecty) - V*t;
						lago[j][k] += aux * exp(-1*aux*aux - (t/10));
					}
				}
			}
		}
		if(rand() <= ((float) Z/100) * RAND_MAX) {
			gotas[numGotas].x = inteiroAleatorio(L);
			gotas[numGotas].y = inteiroAleatorio(H);
			printf("gerou uma gota: %d %d \n", gotas[numGotas].x, gotas[numGotas].y);
			gotas[numGotas].tempo = (float) i * T / NITER;
			numGotas++;
		}
	}
	

	escreveArquivo(lago);
	return 0;
}

float** criaMatriz (int m, int n)
{
    float** M; 
	int i;
	M = malloc (m * sizeof (float *));
	for (i = 0; i < m; ++i)
	   M[i] = calloc (n, sizeof (float));

	return M;
}

/* gera inteiro entre 0 e maximo - 1 */
int inteiroAleatorio (int maximo)
{
    return ( (double) rand() / ((double) RAND_MAX + 1)) * maximo;
}

float calculaDistancia (int x1, int y1, int x2, int y2, float aspectx, float aspecty) {
	float difx = (x2 - x1)*aspectx;
	float dify = (y2 - y1)*aspecty; 
	return sqrt( difx * difx  + dify * dify );
}

void pegaEntrada(int argc, char* argv[]) {
	FILE *entrada;
	if(argc > 1) {
		entrada = fopen(argv[1], "r");
		fscanf(entrada, "(%f,%f) \n (%d,%d) \n %d \n %f \n %f \n %d \n %f \n %d",
			&LARG, &ALT, &L, &H, &T, &V, &EPS, &NITER, &P, &SEED);
		fclose(entrada);

		if(argc > 2)
			NPROCS = atoi(argv[2]);
	}
}

void escreveArquivo(float** lago) {

	int i, j;
	float hmax = 0.0, pmax = 0.0, delta;
	FILE *saida;
	saida = fopen("saida.ppm", "w");
	fprintf(saida, "P3\n");
	fprintf(saida, "%d %d \n255\n", L, H);
	
	for(i = 0; i < H; i++) {
		for(j = 0; j < L; j++) {
			if(lago[i][j] > hmax) {
				hmax = lago[i][j];
			}
			else if(lago[i][j] < pmax) {
				pmax = lago[i][j];
			}
		}
	}

	delta = (hmax > -pmax)? hmax/255 : -pmax/255;  

	for(i = 0; i < H; i++) {
		for(j = 0; j < L; j++) {
			if(lago[i][j] < 0) {
				fprintf(saida, "%d 0 0\n", (int) ceil(-lago[i][j]/delta));
			}
			else {
				fprintf(saida, "0 0 %d\n", (int) ceil(lago[i][j]/delta));
			}
				
		}
	}
	fclose(saida);
}