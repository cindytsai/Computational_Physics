#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Do integration from a to b.
// The points index start from 0 to N, 
// so there are N interval and N+1 points.
double f(double x);
double linear(int N, double a, double b);
double quadratic(int N, double a, double b);
double quartic(int N, double a, double b);

void main(){
	// Print a table like structure, and some info.
	printf("exact solution : 2\n");
	printf("Solution Given By Following Methods=================================================================\n\n");
	printf("  N  \t\t\tlinear\t\t\tquadratic\t\t\tquartic\n");

	// number of intervals
	int N;
	// Find the integral by these methods
	for(int i = 2; i <= 10; i = i+1){
		N = (int)pow(2, i);
		printf("%4d", N);
		printf("\t%1.20e", linear(N, 0.0, M_PI));
		printf("\t%1.20e", quadratic(N, 0.0, M_PI));
		printf("\t%1.20e\n", quartic(N, 0.0, M_PI));
	}

	printf("\nErrors==============================================================================================\n\n");
	printf("  N  \t\t\tlinear\t\t\tquadratic\t\t\tquartic\n");
	// Find the errors for each methods.
	for(int i = 2; i <= 10; i = i+1){
		N = (int)pow(2, i);
		printf("%4d", N);
		printf("\t%1.20e", linear(N, 0.0, M_PI) - 2.0);
		printf("\t%1.20e", quadratic(N, 0.0, M_PI) - 2.0);
		printf("\t%1.20e\n", quartic(N, 0.0, M_PI) - 2.0);
	}
	exit(0);
}

double f(double x){
	return sin(x);
}

double linear(int N, double a, double b){
	// index = 0 1 2 3 4 ... ... N
	// h/2 * ( 1 2 2 2 2 ... ... 1 )
	double sum = 0.0;
	double h;

	h = (b - a) / N;
	sum = sum + f(a) + f(b);
	for(int i = 1; i < N; i = i+1){
		sum = sum + 2 * f(a + h * i);
	}
	sum = (h / 2.0) * sum;

	return sum;
}

double quadratic(int N, double a, double b){
	// index = 0 1 2 3 4 5 6 ... ... N
	// h/3 * ( 1 4 2 4 2 4 2 ... ... 1)
	double sum = 0.0;
	double h;

	h = (b - a) / N;
	sum = sum + f(a) + f(b);
	for(int i = 1; i < N; i = i+1){
		if( i % 2 == 1) sum = sum + 4 * f(a + h * i);
		if( i % 2 == 0) sum = sum + 2 * f(a + h * i);
	}
	sum = (h / 3.0) * sum;

	return sum;
}

double quartic(int N, double a, double b){
	// index   = 0  1  2  3  4  5  6  7  8  ... ... N
	// 2h/45 * ( 7 32 12 32 14 32 12 32 14  ... ... 7)
	double sum = 0.0;
	double h;

	h = (b - a) / N;
	sum = sum + 7 * (f(a) + f(b));
	for(int i = 1; i < N; i = i+1){
		if( i % 4 == 1) sum = sum + 32 * f(a + h * i);
		if (i % 4 == 2) sum = sum + 12 * f(a + h * i);
		if (i % 4 == 3) sum = sum + 32 * f(a + h * i);
		if (i % 4 == 0) sum = sum + 14 * f(a + h * i);
	}
	sum = ((2 * h) / 45) * sum;

	return sum;
}