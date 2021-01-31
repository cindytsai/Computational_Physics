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
	// N = pow(2, n);
	// n = 1, 2 , ..., 16
	// array x --> holds the random number for coordinate. 
	// Store temperary random number in r, y.
	double meanS, sigmaS;
	double meanR, sigmaR;
	double meanM, sigmaM;
	int N;
	double x[10], x_old[10];
	double r;

	// Output to a file
	FILE *output;
	output = fopen("integration.txt", "w");
	fprintf(output, "N SimpleSampling SSerror RejectionMethod RMerror MetropolisAlgorithm MAerror\n");
	
	for(int n = 1; n <= 16; n = n+1){

		// Sample N points
		N = pow(2, n);
		fprintf(output, "%d ", N);

		/*-----Simple Sampling-----*/
		meanS = 0.0;
		sigmaS = 0.0;
		for(int i = 1; i <= N; i = i+1){
			// Get random x coordinates
			for(int j = 0; j < 10; j = j+1){
				x[j] = gsl_rng_uniform(rng);
			}
			meanS = meanS + f(x);
			sigmaS = sigmaS + pow(f(x), 2);
		}
		meanS = meanS / (double) N;
		sigmaS = sqrt(((1.0 / (double) N) * sigmaS + pow(meanS, 2)) / (double) N);
		fprintf(output, "%.5e +- %.5e ", meanS, sigmaS);

		/*-----Rejection Method-----*/
		// Take g(x) = 1, G(x) = y = x, G^-1(y) = x = y
		// so random y ---> x
		meanR = 0.0;
		sigmaR = 0.0;
		for(int i = 1; i <= N; i = i+1){
			while(1){
				// Get random x coordinates
				for(int j = 0; j < 10; j = j+1){
					x[j] = gsl_rng_uniform(rng);
				}
				r = gsl_rng_uniform(rng);
				if (r <= f(x)){
					meanR = meanR + f(x) / w(x);
					sigmaR = sigmaR + pow(f(x) / w(x), 2);
					break;
				}
			}
		}
		meanR = meanR / (double) N;
		sigmaR = sqrt(((1.0 / (double) N) * sigmaR + pow(meanR, 2)) / (double) N);
		fprintf(output, "%.5e +- %.5e ", meanR, sigmaR);

		/*-----Metropolis Algorithm-----*/
		meanM = 0.0;
		sigmaM = 0.0;
		// Get initial x --> x_old
		for(int j = 0; j < 10; j = j+1){
			x_old[j] = gsl_rng_uniform(rng);
		}
		meanM = meanM + f(x_old) / w(x_old);
		sigmaM = sigmaM + pow(f(x_old) / w(x_old), 2);
		// Get the other (N-1) sample points
		for(int i = 2; i <= N; i = i+1){
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
			meanM = meanM + f(x_old) / w(x_old);
			sigmaM = sigmaM + pow(f(x_old) / w(x_old), 2);
		}
		meanM = meanM / (double) N;
		sigmaM = sqrt(((1.0 / (double) N) * sigmaM + pow(meanM, 2)) / (double) N);
		fprintf(output, "%.5e +- %.5e\n", meanM, sigmaM);
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