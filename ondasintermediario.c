#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
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
float calculaDistancia (int x1, int y1, int x2, int y2);
void escreveArquivo(float** lago);
void pegaEntrada(int argc, char* argv[]);
void zeraLago(float** lago, int ymin, int ymax, int xmin, int xmax);
float h(float p, float t);
void calculaQuadradoLago(float** lago, int xi, int xf, int yi, int yf, float x_gota_lago, float y_gota_lago, float t);

typedef struct {
	int x;
	int y;
	float tempo;
} gota;

int H, L, NITER, T, SEED, NPROCS, ALT, LARG;
float aspectx, aspecty;
float P, V, EPS;

float** tem_gota;

int main (int argc, char* argv[]) {

	int i, l, xi, yi, xf, yf, xmin, xmax, ymin, ymax, numGotas = 0;
	float t, auxx, auxy, max, raizdedois, x_gota_lago, y_gota_lago;
	gota * gotas;
	float** lago;

	pegaEntrada(argc, argv);
	srand(SEED);
	aspectx = (float) LARG/L;
	aspecty = (float) ALT/H;
	gotas = malloc (NITER * sizeof(gota));
	lago = criaMatriz (H, L);
	tem_gota = criaMatriz (H, L);
	raizdedois = sqrt(2)/2;
	

	for(i = 0; i < NITER; i++) {
		if(numGotas != 0) {
			xmin = L; xmax = 0; ymin = L; ymax = 0;
			for(l = 0; l < numGotas; l++) {
				t = (float) i * T / NITER - gotas[l].tempo;
				max = (raizdedois + V*t);
				auxy = (int) (max/aspecty);
				auxx = (int) (max/aspectx);

				yi = gotas[l].y - auxy - 1;
				if(yi < 0) yi = 0;	
				yf =  gotas[l].y + auxy + 1;
				if(yf > H) yf = H;
				xi = gotas[l].x - auxx - 1;
				if(xi < 0) xi = 0;	
				xf =  gotas[l].x + auxx + 1;
				if(xf > L) xf = L;

				x_gota_lago = gotas[l].x * aspectx;
				y_gota_lago = gotas[l].y * aspecty;

				while(yi > 0 && h(calculaDistancia(x_gota_lago, yi * aspecty, x_gota_lago, y_gota_lago), t) > EPS) {
					yi--;
				}
				while(yf < H && h(calculaDistancia(x_gota_lago, yf * aspecty, x_gota_lago, y_gota_lago), t) > EPS) {
					yf++;
				}
				while(xi > 0 && h(calculaDistancia(xi * aspectx, y_gota_lago, x_gota_lago, y_gota_lago), t) > EPS) {
					xi--;
				}
				while(xf < L && h(calculaDistancia(xf * aspectx, y_gota_lago, x_gota_lago, y_gota_lago), t) > EPS) {
					xf++;
				}

				if(yi < ymin) ymin = yi;
				if(xi < xmin) xmin = xi; 
				if(yf > ymax) ymax = yf;
				if(xf > xmax) xmax = xf;
				//printf("%d %d %d %d \n", xi, yi, xf, yf);

				calculaQuadradoLago(lago, xi, xf, yi, yf, x_gota_lago, y_gota_lago, t);
			}
			if(i < NITER -1)
			zeraLago(lago, ymin, ymax, xmin, xmax);
		}
		if(rand() <= ((float) P/100) * RAND_MAX) {
			gotas[numGotas].x = inteiroAleatorio(L);
			gotas[numGotas].y = inteiroAleatorio(H);

			tem_gota[(int) (gotas[numGotas].y)][(int) (gotas[numGotas].x)] = 1;

			printf("gerou uma gota de numero %d: %d %d \n", numGotas, gotas[numGotas].x, gotas[numGotas].y);
			gotas[numGotas].tempo = (float) i * T / NITER;
			numGotas++;

		}
	}
	
	// for(j = 0; j < H; j++) {
	// 	for(k = 0; k < L; k++) {
	// 		printf("%f ", lago[j][k]);
	// 	}
	// 	printf("\n");
	// }
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

float calculaDistancia (int x1, int y1, int x2, int y2) {
	float difx = (x2 - x1);
	float dify = (y2 - y1); 
	return sqrt( difx * difx  + dify * dify );
}

void pegaEntrada(int argc, char* argv[]) {
	FILE *entrada;
	if(argc > 1) {
		entrada = fopen(argv[1], "r");
		fscanf(entrada, "(%d,%d) \n (%d,%d) \n %d \n %f \n %f \n %d \n %f \n %d",
			&LARG, &ALT, &L, &H, &T, &V, &EPS, &NITER, &P, &SEED);
		fclose(entrada);

		if(argc > 2)
			NPROCS = atoi(argv[2]);
	}
	else
	{
		printf("Uso: %s <entrada>\n", argv[0]);
		exit(EXIT_FAILURE);
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

	printf("hmax: %f pmax: %f \n", hmax, pmax); 
	delta = (hmax > -pmax)? hmax/255 : -pmax/255;  
	if(delta == 0)
		delta = 1;
	for(i = 0; i < H; i++) {
		for(j = 0; j < L; j++) {
			if(tem_gota[i][j] == 1)
				fprintf(saida, "0 255 0\n");
			else if(tem_gota[i][j] == 2)
				fprintf(saida, "255 255 255\n");
			else
			{
				if(lago[i][j] < 0) {
					fprintf(saida, "%d 0 0\n", (int) ceil(-lago[i][j]/delta));
				}
				else {
					fprintf(saida, "0 0 %d\n", (int) ceil(lago[i][j]/delta));
				}
			}			
		}
	}
	fclose(saida);
}

void zeraLago(float** lago, int ymin, int ymax, int xmin, int xmax) {
	int i, j;
	for(i = ymin; i < ymax; i++)
		for(j = xmin; j < xmax; j++)
			lago[i][j] = 0;
}

float h(float ro, float t) {
	float aux = ro - V*t;
	return aux * exp(-1*aux*aux - (t/10));
}

void calculaQuadradoLago(float** lago, int xi, int xf, int yi, int yf, float x_gota_lago, float y_gota_lago, float t) {
	int j, k;
	float altura;
	#pragma omp parallel for private(altura, k) 
	for(j = yi; j < yf; j++) {
		for(k = xi; k < xf; k++) {
			altura = h(calculaDistancia(k * aspectx, j * aspecty, x_gota_lago, y_gota_lago), t);
			
			if(fabs(altura) >= EPS) {
				lago[j][k] += altura;
				
			}
		}
	}
}