/*
Use Richardson Extrapolation to find second derivatives of f(x) = x * exp(x)
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);
double f2prime(double x);
double threePoint(double x, double h);

void main(){

	// Settings: h is initial step , and calculate with N interations
	double h = 0.4;
	double x = 2.0;
	int N = 4;

	// Print out a table like structure, and exact solution.
	printf("f(x) = x * exp(x)\n");
	printf("f''(x) = 2 * exp(x) + x * exp(x)\n");
	printf("f''(2) = %1.10e\n", f2prime(x));
	printf("============================================================================\n");
	printf("Using Richardson Extrapolation Method:\n");
	printf(" h\t\tD1\t\tD2\t\tD3\t\tD4\n");

	// Store richardson extra. in D matrix.
	double D[N][N];
	
	// Calculating D, and print it.
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){

			if(i < j) break;

			else if( j == 0){
				D[i][j] = threePoint(x, h);
				printf("%1.2f  ", h);
				printf("%1.10e  ", D[i][j]);
			}

			else{
				D[i][j] = (pow(4.0, j) * D[i][j-1] - D[i-1][j-1]) / (pow(4.0, j) - 1);
				printf("%1.10e  ", D[i][j]);
			}

		}
		printf("\n");
		h = h / 2.0;
	}

}

double f(double x){
	return x * exp(x);
}

double f2prime(double x){
	return 2.0 * exp(x) + x * exp(x);
}

double threePoint(double x, double h){
	return (f(x+h) - 2.0 * f(x) + f(x-h)) / pow(h, 2);
}