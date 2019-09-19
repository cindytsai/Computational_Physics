#include<stdio.h>
#include<stdlib.h>

void main(){
		double x; 
		double eps = 1.0;
		int i, N = 100;
		for (i = 1; i < N; i = i + 1){
				eps = eps / 2.0;
				x = 1.0 + eps;
				printf("i = %d, x = %25.20e, eps = %25.20e \n", i, x, eps);

				if (x == 1.0) {
						exit(0);
				}
		}
}
