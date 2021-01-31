#include<stdio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<string.h>
#include <stdlib.h>

// Global
// J   -> Coupling strength
// B   -> External magnetic field
// T   -> Temperature
// kB  -> Boltzmann constant
double J = 1.0;
double B = 0.0;
double kB = 1.0;

// Functions
// Output array as a file
void outputFile(double *array_e, double *array_m, char *file_name, int array_length);

// Algorithms
void metropolis(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name);
void heat_bath(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name);
void single_cluster(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name);


int main(void){
	// Settings and Notations 
	// nt  -> numbers of sweeps for thermalization
	// nm  -> numbers of measurements, and make sure it is 2^n
	// im  -> interval between successive measurements
	// seed    -> seed for random number generator
	// *rng    -> random number generator
	// 
	// Tc_top / Tc_bottom -> Tc should have be inside this boundary
	// T                  -> Store temperary temperature
	// ave_m              -> Store <m> under that temperature
	// jump               -> Record if <m> jumps over 0.47 < 0.53
	//                    -> If true => jump = 1;
	// input              -> Read file during process
	int nt = 20000;
	int nm = 10000;
	int im = 1;
	int seed = 4357;
	gsl_rng *rng;
	double Tc_bottom = 0.000, Tc_top = 2.701;
	double T = Tc_top; 
	double ave_m = 0.0;
	int jump = 0;
	char file_name[200], temp_file[200], temp_str[50];
	FILE *input, *output;
	double temp, m;
	int file_operation;


	// Initializing
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);

	output = fopen("L-Tc.txt", "a");
	fprintf(output, "L Tc\n");
	fclose(output);
	

	for(int L = 10; L <= 200; L = L+10){
		// ns      -> total # of spin
		// spin    -> spin array for the 2D lattice
		int ns = L * L;
		int spin[ns];
		// Create a spin lattice
		for(int i = 0; i < ns; i = i+1){
			if(gsl_rng_uniform(rng) > 0.5) spin[i] = 1;
			else spin[i] = -1;
		}

		// Initialized
		if (jump == 1){
			T = T - 0.0005 + 0.01 ;
		}
		else{
			T = T + 0.01;
		}
		jump = 0;
		ave_m = 0.0;

		printf("L : %d\n", L);

		// Find Tc in range T_bottom~T_top
		//      L              200     10
		// Break while loop if <m> suddenly jumps very large
		while((fabs(ave_m - 0.5) > 0.03) && (T >= Tc_bottom) && (jump == 0)){

			T = T - 0.001;
			ave_m = 0.0;
			file_operation = remove("temp.txt");

			printf("T : %.3e\n", T);

			// Set temp file name
			strcpy(temp_file, "temp.txt");
			// Call Algorithm, do thermalization, measurements, and output 
			// metropolis(spin, L, T, nt, nm, im, rng, temp_file);
			// heat_bath(spin, L, T, nt, nm, im, rng, temp_file);
			single_cluster(spin, L, T, nt, nm, im, rng, temp_file);
			
			// Calculate <m> in file temp.txt
			input = fopen(temp_file, "r");
			for(int i = 0; i < nm; i = i+1) {
      			fscanf(input,"%lf %lf", &temp, &m);
      			ave_m = ave_m + m;
    		}
    		fclose(input);
    		ave_m = ave_m / nm;

    		if(ave_m >= 0.53){
    			jump = 1;
    			T = T + 0.0005;
    		}

    		printf("<m> = %.3e \n", ave_m);
		}

		// Now T = Tc, output to file L-Tc.txt
		output = fopen("L-Tc.txt", "a");
		fprintf(output, "%d %.4e\n", L, T);
		fclose(output);
		printf("L: %d, Tc = %.4e\n", L, T);

		// Set file name
		strcpy(file_name, "");
		strcat(file_name, "L");
		sprintf(temp_str, "%d", L);
		strcat(file_name, temp_str);
		strcat(file_name, ".txt");

		// Rename temp.txt to LXXX.txt
		file_operation = rename(temp_file, file_name);
		if(file_operation == 0){
			printf("L%d Rename Successful !\n", L);
		}
		else{
			printf("L%d Rename Failed !\n", L);
		}		
	}
	

	return 0;
}

void outputFile(double *array_e, double *array_m, char *file_name, int array_length){
	// Create txt file file_name
	FILE *output;
	output = fopen(file_name, "w");
	// Print results
	for(int i = 0; i < array_length; i = i+1){
		fprintf(output, "%1.5e %1.5e\n", array_e[i], array_m[i]);
	}
	// Close file
	fclose(output);
}

void metropolis(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name){
	// Settings
	// sweeps    -> total # of sweeps
	// swp       -> record sweep numbers
	// fw/bw     -> forward and backward index
	// k         -> spin index
	// k1~k4     -> neighboring spin index [right, top, left, bottom]
	// old/new   -> old/new spin sigma
	// spins     -> sum of neighboring spin sigma
	// dE        -> delta E = E(new) - E(old)
	// mag       -> total magnetization
	// energy    -> total energy
	int sweeps = nt + nm * im;
	int swp;
	int fw[N], bw[N];
	int k, k1, k2, k3, k4;
	int old, new;
	int spins;
	double dE;
	double mag, energy;
	// Store individual measurements
	// e  -> energy density
	// m  -> magnetization density
	// E2 -> energy square
	// M2 -> magnetization square
	int store_index = 0;
	int array_swp[nm];
	double array_e[nm], array_m[nm], array_E2[nm], array_M2[nm];

	for(int i = 0; i < N; i = i+1) {
		fw[i] = (i + 1) % N;
		bw[i] = (i - 1 + N) % N;
	}

	// Sweep through lattice, and measure
	for(swp = 1; swp <= sweeps; swp = swp+1){

		// Thermalization for every iteration
		for(int j = 0; j < N; j = j+1){
			for(int i = 0; i < N; i = i+1){
				k = i + j * N;
				old = spin[k];
				new = -old;
				k1 = fw[i] + j * N;
				k2 = i + fw[j] * N;
				k3 = bw[i] + j * N;
				k4 = i + bw[j] * N;
				spins = spin[k1] + spin[k2] + spin[k3] + spin[k4];
				dE = -(new - old) * (J * spins + B);

				if((dE <= 0.0) || (gsl_rng_uniform(rng) <= exp(-dE / (kB * T)))){
					spin[k] = new;
				}
			}
		}

		// Measure
		if((swp > nt) && (swp % im == 0)){
			mag = 0.0;
			energy = 0.0;

			for(int j = 0; j < N; j = j+1){
				for(int i = 0; i < N; i = i+1){
					k = i + j * N;
					k1 = fw[i] + j * N;
					k2 = i + fw[j] * N;

					mag = mag + spin[k];
					energy = energy - J * spin[k] * (spin[k1] + spin[k2]);
				}
			}
			energy = energy - B * mag;

			// Store measurement to array
			array_swp[store_index] = swp;
			array_e[store_index] = energy / (N * N);
			array_m[store_index] = fabs(mag) / (N * N);
			array_E2[store_index] = pow(energy, 2);
			array_M2[store_index] = pow(mag, 2);
			store_index = store_index + 1;
		}
	}

	// Output result to file
	outputFile(array_e, array_m, file_name, nm);
}

void heat_bath(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name){
	// Settings
	// sweeps    -> total # of sweeps
	// swp       -> record sweep numbers
	// fw/bw     -> forward and backward index
	// k         -> spin index
	// k1~k4     -> neighboring spin index [right, top, left, bottom]
	// spins     -> sum of neighboring spin sigma
	// dE        -> delta E = E(new) - E(old)
	// mag       -> total magnetization
	// energy    -> total energy
	int sweeps = nt + nm * im;
	int swp;
	int fw[N], bw[N];
	int k, k1, k2, k3, k4;
	int spins;
	double dE;
	double mag, energy;
	// Store individual measurements
	// e  -> energy density
	// m  -> magnetization density
	// E2 -> energy square
	// M2 -> magnetization square
	int store_index = 0;
	int array_swp[nm];
	double array_e[nm], array_m[nm], array_E2[nm], array_M2[nm];

	for(int i = 0; i < N; i = i+1) {
		fw[i] = (i + 1) % N;
		bw[i] = (i - 1 + N) % N;
	}

	// Sweep through lattice and measure
	for(swp = 1; swp <= sweeps; swp = swp+1){

		// Thermalization for every iteration
		for(int j = 0; j < N; j = j+1){
			for(int i = 0; i < N; i = i+1){
				k = i + j * N;
				k1 = fw[i] + j * N;
				k2 = i + fw[j] * N;
				k3 = bw[i] + j * N;
				k4 = i + bw[j] * N;
				spins = spin[k1] + spin[k2] + spin[k3] + spin[k4];
				dE = -2 * (J * spins + B);

				if(gsl_rng_uniform(rng) < 1.0 / (1.0 + exp(dE / (kB * T)))) spin[k] = 1;
				else spin[k] = -1;
			}
		}

		// Measure
		if((swp > nt) && (swp % im == 0)){
			mag = 0.0;
			energy = 0.0;

			for(int j = 0; j < N; j = j+1){
				for(int i = 0; i < N; i = i+1){
					k = i + j * N;
					k1 = fw[i] + j * N;
					k2 = i + fw[j] * N;

					mag = mag + spin[k];
					energy = energy - J * spin[k] * (spin[k1] + spin[k2]);
				}
			}
			energy = energy - B * mag;

			// Store measurement to array
			array_swp[store_index] = swp;
			array_e[store_index] = energy / (N * N);
			array_m[store_index] = fabs(mag) / (N * N);
			array_E2[store_index] = pow(energy, 2);
			array_M2[store_index] = pow(mag, 2);
			store_index = store_index + 1;
		}
	}


	// Output result to file
	outputFile(array_e, array_m, file_name, nm);
}

void single_cluster(int *spin, int N, double T, int nt, int nm, int im, gsl_rng *rng, char *file_name){
	// Settings
	// sweeps    -> total # of sweeps
	// swp       -> record sweep numbers
	// fw/bw     -> forward and backward index
	// i1, j1    -> temp coordinate
	// old/new   -> old/new spin sigma
	// k         -> spin index
	// k1, k2    -> neighboring spin index [right, top]
	// kk[4]     -> neighboring spin index [right, top, left, bottom]
	// spins     -> sum of neighboring spin sigma
	// pt        -> pointer for cluster
	// cluster   -> cluster of spin index per one iteration
	// prob      -> probability of forming bond
	// dE        -> delta E = E(new) - E(old)
	// mag       -> total magnetization
	// energy    -> total energy
	int sweeps = nt + nm * im;
	int swp;
	int fw[N], bw[N];
	int i1, j1;
	int old, new;
	int k, k1, k2;
	int kk[4];
	int spins;
	int pt;
	int cluster[N * N];
	double prob = 1.0 - exp(-2 * J / (kB * T));
	double dE;
	double mag, energy;
	// Store individual measurements
	// e  -> energy density
	// m  -> magnetization density
	// E2 -> energy square
	// M2 -> magnetization square
	int store_index = 0;
	int array_swp[nm];
	double array_e[nm], array_m[nm], array_E2[nm], array_M2[nm];

	for(int i = 0; i < N; i = i+1) {
		fw[i] = (i + 1) % N;
		bw[i] = (i - 1 + N) % N;
	}

	// Sweep through lattice and measure
	for(swp = 1; swp <= sweeps; swp = swp+1){

		// Thermalization for every iteration
		for(int i = 0; i < N / 4; i = i+1){
			
			// Initialize cluster
			for(int j = 0; j < N * N; j = j+1){
				cluster[j] = 0;
			}
			
			// Pick a random site, and flip the spin
			i1 = ((int)(N * gsl_rng_uniform(rng))) % N;
			j1 = ((int)(N * gsl_rng_uniform(rng))) % N;
			// i1 = N * gsl_rng_uniform(rng);
			// j1 = N * gsl_rng_uniform(rng);
			// i1 = i1 % N;
			// j1 = j1 % N;
			k = i1 + j1 * N;
			old = spin[k];
			new = -old;
			spin[k] = new;

			// Push the site to cluster, and loop through cluster until it's not empty
			pt = 0;
			cluster[pt] = k;
			while((pt >= 0) && (cluster[pt] >= 0)){
				k = cluster[pt];
				cluster[pt] = 0;
				pt = pt - 1;
				i1 = k % N;
				j1 = (k - i1) / N;
				kk[0] = fw[i1] + j1 * N;
				kk[1] = i1 + fw[j1] * N;
				kk[2] = bw[i1] + j1 * N;
				kk[3] = i1 + bw[j1] * N;

				for(int j = 0; j < 4; j = j+1){
					
					// Accept the neighboring spin, and add it to cluster
					if((spin[kk[j]] == old) && (gsl_rng_uniform(rng) < prob)){
						spin[kk[j]] = new;
						pt = pt + 1;

						// In case of error
						if((pt >= N * N) || (pt < 0)){
							printf("pt = %d\n", pt);
							exit(1);
						}

						cluster[pt] = kk[j];
					}
				}
			}
		}

		// Measure
		if((swp > nt) && (swp % im == 0)){
			mag = 0.0;
			energy = 0.0;

			for(int j = 0; j < N; j = j+1){
				for(int i = 0; i < N; i = i+1){
					k = i + j * N;
					k1 = fw[i] + j * N;
					k2 = i + fw[j] * N;

					mag = mag + spin[k];
					energy = energy - J * spin[k] * (spin[k1] + spin[k2]);
				}
			}
			energy = energy - B * mag;

			// Store measurement to array, and store magnetization as positive
			array_swp[store_index] = swp;
			array_e[store_index] = energy / (N * N);
			array_m[store_index] = fabs(mag) / (N * N);
			array_E2[store_index] = pow(energy, 2);
			array_M2[store_index] = pow(mag, 2);
			store_index = store_index + 1;
		}
	}


	// Output result to file
	outputFile(array_e, array_m, file_name, nm);
}
