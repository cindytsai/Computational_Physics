#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<math.h>
#include<string.h>

double f(double *x);
double w(double *x);

int main(void){
	// Setting up random number generator
	gsl_rng *rng;
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 101);

	// Settings
	// Sampling N points
	// array x --> holds the random number for coordinate. 
	// Store temperary random number in r, y.
	int N = pow(2, 16);
	double x[10], x_old[10];
	double r;

	// Output to a file
	FILE *output;
	output = fopen("data.txt", "w");

	for(int n = 1; n <= N; n = n+1){

		/*-----Simple Sampling-----*/
		// Get random x coordinates
		for(int j = 0; j < 10; j = j+1){
			x[j] = gsl_rng_uniform(rng);
		}
		fprintf(output, "%.5e ", f(x));



		/*-----Metropolis Algorithm-----*/
		// Get initial x --> x_old
		if(n == 1){
			for(int j = 0; j < 10; j = j+1){
				x_old[j] = gsl_rng_uniform(rng);
			}
			fprintf(output, "%.5e\n", f(x_old) / w(x_old));
		}
		else{
			// Get the other (N-1) sample points
			// Get new x --> x
			for(int j = 0; j < 10; j = j+1){
				x[j] = gsl_rng_uniform(rng);
			}
			
			// Check acceptance
			if(w(x) >= w(x_old)){ 
				// Accept x, and to avoid overflow
				memcpy(x_old, x, sizeof(x_old));
			}
			else{
				r = gsl_rng_uniform(rng);
				if(r < (w(x) / w(x_old))){
					// Accept x, and to avoid overflow
					memcpy(x_old, x, sizeof(x_old));
				}
			}
			fprintf(output, "%.5e\n", f(x_old) / w(x_old));
		}
	}

	fclose(output);
	return 0;
}

double f(double *x){
	double result = 1.0;
	for(int i = 0; i <= 9; i = i+1){
		result = result + pow(x[i], 2);
	}
	result = 1.0 / result;
	return result;
}

double w(double *x){
	double result = 0.0;
	for(int i = 0; i <= 9; i= i+1){
		result = result + x[i];
	}
	result = 2.0 - (1.0/5.0) * result;
	return result;
}