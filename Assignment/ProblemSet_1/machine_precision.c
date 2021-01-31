/*
Determine the machine precision of your computer,
with different c data type stored respectively.

float  --> single precision 
double --> double precision
 */

# include<stdio.h>
# include<stdlib.h>

void main(){
	
	float single_prec;
	double double_prec;

	double eps;
	int N;

	// Find epsilon, for single_prec (float)
	eps = 1.0;
	N = 100;

	for(int i = 1; i < N; i = i+1){
		eps = eps / 2.0;
		single_prec = 1.0 + eps;
		if (single_prec == 1.0){
			printf("single precision (float),\nmachine precision epsilon : %5.10e \t (2^-%d) \n", eps, i);
			break;
		}
	}

	// Find epsilon, for double_prec (double)
	eps = 1.0;
	N = 100;

	for (int i = 1; i < N; i = i+1){
		eps = eps / 2.0;
		double_prec = 1.0 + eps;
		if (double_prec == 1.0){
			printf("double precision (double),\nmachine precision epsilon : %5.10e \t (2^-%d) \n", eps, i);
			exit(0);
		}
	}
}