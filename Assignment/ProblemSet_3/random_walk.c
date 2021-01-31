#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<math.h>

int main(void){

	// Setting up random number generator
	gsl_rng *rng;
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 1000);

	// Output to a file
	FILE *output;
	output = fopen("output.txt", "w");
	fprintf(output, "N\t\td\n");

	// Settings
	// r --> holds the random number
	// from 1 ~ N total steps, and each total steps simulate total of M times.
	// coor for recording the coordinates, origin starts at (0, 0)
	// sum up M times of the result distance.
	unsigned long int r;
	int N = 10000;
	int M = 100;
	int coor[2];
	double sum_d;

	for(int i = 1; i <= N; i = i+1){

		sum_d = 0.0;

		for(int j = 1; j <= M; j = j+1){

			coor[0] = 0;
			coor[1] = 0;

			// Move N steps
			for(int k = 1; k <= i; k = k+1){
				r = gsl_rng_get(rng);
				if(r % 4 == 0) coor[0] = coor[0] + 1;
				if(r % 4 == 1) coor[1] = coor[1] + 1;
				if(r % 4 == 2) coor[0] = coor[0] - 1;
				if(r % 4 == 3) coor[1] = coor[1] - 1;
			}
			sum_d = sum_d + sqrt(pow(coor[0], 2) + pow(coor[1], 2));
		}

		fprintf(output, "%d\t\t%e\n", i, sum_d / M);
	}

	fclose(output);
	return 0;
}