#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(char c, double w);
double linear(char c, int N, double v);
double rombergIntergration(char c, double v);
double I_I0(double v);

void main(){
	// Create output file
	FILE *output;
	output = fopen("output.txt", "w");

	// Print a table like structure
	fprintf(output, "   v   \tI/I0\n");

	// Setting
	// Find I/I0 for v range 0 ~ 100
	// with v_N interval, v_N + 1 points
	double v_max = 100.0;
	int v_N = 10000;
	double steps = v_max / v_N;

	// Calculate I/I0, with different v
	for(int i = 0; i <= v_N; i = i+1){
		fprintf(output, "%3.2f", i * steps);
		fprintf(output, "   %e\n", I_I0(i * steps));
	}

	fclose(output);
	exit(0);
}

double f(char c, double w){
	// There are two cases
	if(c == 'C') return cos(M_PI * pow(w, 2) * 0.5);
	if(c == 'S') return sin(M_PI * pow(w, 2) * 0.5);
}

double linear(char c, int N, double v){
	// index = 0 1 2 3 4 ... ... N
	// h/2 * ( 1 2 2 2 2 ... ... 1 )
	double sum = 0.0;
	double h;

	h = v / N;
	sum = sum + f(c, 0) + f(c, v);
	for(int i = 1; i < N; i = i+1){
		sum = sum + 2 * f(c, h * i);
	}
	sum = (h / 2.0) * sum;

	return sum;
}

double rombergIntergration(char c, double v){
	// Do the intergration from 0 to v.
	// with interval n = 2^10, 2^11 2^12 2^13
	// with k interation.
	// Store it in double matrix T
	int k = 4;
	int m = 10;
	double T[k][k];

	for(int i = 0; i < k; i = i+1){
		for(int j = 0; j < k; j = j+1){

			if(i < j) break;
			else if (j == 0){
				T[i][j] = linear(c, pow(2, m+i), v);
			}
			else {
				T[i][j] = (pow(4, j) * T[i][j-1] - T[i-1][j-1]) / (pow(4, j) - 1);
			}
		}
	}

	return T[k-1][k-1];
}

double I_I0(double v){

	double intensity = 0.0;

	intensity = intensity + pow(rombergIntergration('C', v) + 0.5, 2);
	intensity = intensity + pow(rombergIntergration('S', v) + 0.5, 2);
	intensity = 0.5 * intensity;

	return intensity;
}