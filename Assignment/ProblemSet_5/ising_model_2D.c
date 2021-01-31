#include<stdio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<string.h>
#include <stdlib.h>

// Global
// N   -> N sites per row and column
// J   -> Coupling strength
// B   -> External magnetic field
// T   -> Temperature
// kB  -> Boltzmann constant
int N = 100;
double J = 1.0;
double B = 0.0;
double T = 2.26;
double kB = 1.0;

// Functions
// Output array as a file
void outputFile(int *array_swp, double *array_e, double *array_m, double *array_E2, double *array_M2, char *file_name, int array_length);
// Binning with Jackknife method
void getBinningJackknife(double *A, int bmax2, int N, char *file_name);
// Algorithms
void metropolis(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng);
void heat_bath(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng);
void single_cluster(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng);


int main(void){
	// Settings and Notations 
	// nt  -> numbers of sweeps for thermalization
	// nm  -> numbers of measurements, and make sure it is 2^n
	// im  -> interval between successive measurements
	// seed    -> seed for random number generator
	// bmax2   -> for jackknife binning method 2^1~2^bmax2
	int bmax2 = 10;
	int nt = 1000;
	int nm = pow(2, bmax2);
	int im = 10;
	int seed = 100;
	// ns      -> total # of spin
	// spin    -> spin array for the 2D lattice
	// *rng    -> random number generator
	int ns = N * N;
	int spin[ns];
	gsl_rng *rng;

	// Initializing
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	for(int i = 0; i < ns; i = i+1){
		if(gsl_rng_uniform(rng) > 0.5) spin[i] = 1;
		else spin[i] = -1;
	}

	// Call Algorithm, do thermalization, measurements, binning with jackknife, and output
	metropolis(spin, nt, nm, im, bmax2, rng);
	heat_bath(spin, nt, nm, im, bmax2, rng);
	single_cluster(spin, nt, nm, im, bmax2, rng);

	return 0;
}

void outputFile(int *array_swp, double *array_e, double *array_m, double *array_E2, double *array_M2, char *file_name, int array_length){
	// Create txt file file_name
	FILE *output;
	output = fopen(file_name, "w");
	// Print results
	fprintf(output, "sweep    e    m    E2    M2\n"); 
	for(int i = 0; i < array_length; i = i+1){
		fprintf(output, "%d    %1.5e    %1.5e    %1.5e    %1.5e\n", array_swp[i], array_e[i], array_m[i], array_E2[i], array_M2[i]);
	}
	// Close file
	fclose(output);
}

void getBinningJackknife(double *A, int bmax2, int N, char *file_name){
    // General Settings
    // numbers of blocks m from : 2^1 ~ 2^bmax2
    // numbers of data per block --> nb
    // block average --> B[i]
    // # of array A  --> N
    int m, nb;
    double *B, *Bjk;
    double Bjk_ave, errorJK;

    // Output as file with file name file_name
    FILE *output;
    output = fopen(file_name, "w");

    fprintf(output, "m    nb    Bjk_ave    errorJK\n");

    for(int i = bmax2; i >= 1; i = i-1){
        
        m = pow(2, i);
        nb = N / m;
        B = (double*) malloc(m * sizeof(double));
        Bjk = (double*) malloc(m * sizeof(double));
        Bjk_ave = 0.0;
        errorJK = 0.0;

        // Calculate B_i
        for(int j = 0; j < m; j = j+1){
            B[j] = 0.0;
        }
        for(int j = 0; j < N; j = j+1){
            B[j/nb] = B[j/nb] + (A[j] / (double)nb);
        }

        // Calculate Bjk_s
        for(int s = 0; s < m; s = s+1){
            Bjk[s] = 0.0;
            for(int j = 0; j < m; j = j+1){
                if(j != s){
                    Bjk[s] = Bjk[s] + B[j];
                }
            }
            Bjk[s] = Bjk[s] / (double)(m-1);
        }

        // Calculate Bjk_ave
        for(int j = 0; j < m; j = j+1){
            Bjk_ave = Bjk_ave + Bjk[j];
        }
        Bjk_ave = Bjk_ave / (double)m;

        // Calculate errorJK
        for(int j = 0; j < m; j = j+1){
            errorJK = errorJK + pow((Bjk[j] - Bjk_ave), 2);
        }
        errorJK = errorJK * ((double)(m - 1) / (double)m);
        errorJK = sqrt(errorJK);

        // Print out the result
        fprintf(output, "%d    %d    %1.5e    %1.5e\n", m, nb, Bjk_ave, errorJK);

    }
    // Close file
    fclose(output);
}

void metropolis(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng){
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
	// Output file name
	char file_name[200];

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
	strcpy(file_name, "metropolis_thermalized.txt");
	outputFile(array_swp, array_e, array_m, array_E2, array_M2, file_name, nm);
	// Binning with Jackknife method
	strcpy(file_name, "metropolis_thermalized_e_errorJK.txt");
	getBinningJackknife(array_e, bmax2, nm, file_name);
	strcpy(file_name, "metropolis_thermalized_m_errorJK.txt");
	getBinningJackknife(array_m, bmax2, nm, file_name);
	strcpy(file_name, "metropolis_thermalized_E2_errorJK.txt");
	getBinningJackknife(array_E2, bmax2, nm, file_name);
	strcpy(file_name, "metropolis_thermalized_M2_errorJK.txt");
	getBinningJackknife(array_M2, bmax2, nm, file_name);
}

void heat_bath(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng){
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
	// Output file name
	char file_name[200];

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
	strcpy(file_name, "heatbath_thermalized.txt");
	outputFile(array_swp, array_e, array_m, array_E2, array_M2, file_name, nm);
	// Binning with Jackknife method
	strcpy(file_name, "heatbath_thermalized_e_errorJK.txt");
	getBinningJackknife(array_e, bmax2, nm, file_name);
	strcpy(file_name, "heatbath_thermalized_m_errorJK.txt");
	getBinningJackknife(array_m, bmax2, nm, file_name);
	strcpy(file_name, "heatbath_thermalized_E2_errorJK.txt");
	getBinningJackknife(array_E2, bmax2, nm, file_name);
	strcpy(file_name, "heatbath_thermalized_M2_errorJK.txt");
	getBinningJackknife(array_M2, bmax2, nm, file_name);
}

void single_cluster(int *spin, int nt, int nm, int im, int bmax2, gsl_rng *rng){
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
	// Output file name
	char file_name[200];

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
	strcpy(file_name, "singlecluster_thermalized.txt");
	outputFile(array_swp, array_e, array_m, array_E2, array_M2, file_name, nm);
	// Binning with Jackknife method
	strcpy(file_name, "singlecluster_thermalized_e_errorJK.txt");
	getBinningJackknife(array_e, bmax2, nm, file_name);
	strcpy(file_name, "singlecluster_thermalized_m_errorJK.txt");
	getBinningJackknife(array_m, bmax2, nm, file_name);
	strcpy(file_name, "singlecluster_thermalized_E2_errorJK.txt");
	getBinningJackknife(array_E2, bmax2, nm, file_name);
	strcpy(file_name, "singlecluster_thermalized_M2_errorJK.txt");
	getBinningJackknife(array_M2, bmax2, nm, file_name);
}
