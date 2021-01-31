#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void getAutocorrelation(double *A, int tmax, int N);
struct Binning getBinningMethod(double *A, int bmax2, int N);

struct Binning{
    double average;
    double error;
    // Split data into m blocks
    int m;
};

int main(int argc, char *argv[]){
	// Settings
	// N successive measurements
	// A[1~N] for storing data --> Simple Sampling
	// B[1~N] for storing data --> Metropolis Algorithm
    // 
	// t for computing the integrated autocorrelation time, from t : 0 ~ tmax
    // bmax2 for max of the numbers of blocks, 2^bmax2
	int N = pow(2, 10);
	double *A, *B;
    int tmax = 10000;
    int bmax2 = 10;
    struct Binning result;

    A = (double*) malloc((N + 1) * sizeof(double));
    B = (double*) malloc((N + 1) * sizeof(double));

	// Read from data.txt
	FILE *infile;
	infile = fopen(argv[1], "r");
	for(int i = 1; i <= N; i = i+1) {
      fscanf(infile,"%lf %lf", &A[i], &B[i]);
    }
    fclose(infile);


    // Print out the result from binning 
    if(atoi(argv[2]) == 1){
        result = getBinningMethod(A, bmax2, N);
    }
    if(atoi(argv[2]) == 2){
        result = getBinningMethod(B, bmax2, N);
    }
    printf("%lf %lf %d\n", result.average, result.error, result.m);

    return 0;
}

// void getAutocorrelation(double *A, int tmax, int N){
//     // General settings
//     double A_ave = 0.0, A2_ave = 0.0, errA = 0.0;
//     double *A_cor;
//     double tau;
//     double C;
//     int num;
//     int sat;

//     A_cor = (double*) malloc((tmax + 1) * sizeof(double));

//     // Calculate autocorrelation
//     for(int i = 1; i <= N; i = i+1){
//         A_ave = A_ave + A[i];
//         A2_ave = A2_ave + pow(A[i], 2);
//     }
//     A_ave = A_ave / (double)N;
//     A2_ave = A2_ave / (double)N;
//     errA = sqrt((A2_ave - pow(A_ave, 2)) / (double)(N-1));

//     for(int j = 0; j <= tmax; j = j+1){
//         A_cor[j] = 0.0;
//         num = 0;
//         for(int i = 1; i <= N-j; i = i+1){
//             A_cor[j] = A_cor[j] + A[i] * A[i+j];
//             num = num + 1;
//         }
//         A_cor[j] = (A_cor[j] / (double)num) - pow(A_ave, 2);
//     }

//     // Find integrated autocorrelation time
//     // See if it reach saturation, reach --> sat = 1
//     tau = 0.5;
//     sat = 0;
//     for(int t = 1; t <= tmax; t = t+1){
//         C = A_cor[t] / A_cor[0];
//         if(C < 0.0){
//             printf("t = %d, ", t);
//             sat = 1;
//             break;
//         }
//         tau = tau + C;
//     }

//     if(sat != 1){
//         printf("Not saturate !\n");
//     }
//     else{
//         printf("tau_int = %.5e, sqrt(2 * tau_int) = %.5e\n", tau, sqrt(2.0 * tau));
//         printf("<A> = %.5e +- %.5e\n\n", A_ave, errA * sqrt(2.0 * tau));
//     }
// }

struct Binning getBinningMethod(double *A, int bmax2, int N){
    // General Settings
    // numbers of blocks m from : 2^1 ~ 2^bmax2
    // numbers of data per block --> nb
    // block average --> B[i]
    // 
    // Store errors as array errors, with different binning size
    int m, nb;
    double B_ave, B2_ave, errB;
    double *B;
    double *errors;
    struct Binning result = {.average = 0.0, .error = 0.0, .m = pow(2, bmax2)};

    errors = (double*) malloc((bmax2 + 1) * sizeof(double));

    // printf("  m      nb       1/nb        deltaA\n");

    for(int i = bmax2; i >= 1; i = i-1){
        
        m = pow(2, i);
        nb = N / m;
        B = (double*) malloc(m * sizeof(double));

        for(int j = 0; j < m; j = j+1){
            B[j] = 0.0;
        }
        for(int j = 0; j < N; j = j+1){
            B[j/nb] = B[j/nb] + (A[j+1] / (double)nb);
        }

        B_ave = 0.0;
        B2_ave = 0.0;
        for(int j = 0; j < m; j = j+1){
            B_ave = B_ave + B[j];
            B2_ave = B2_ave + pow(B[j], 2);
        }
        B_ave = B_ave / m;
        B2_ave = B2_ave / m;
        errB = sqrt((B2_ave - pow(B_ave, 2)) / (double)(m - 1));

        errors[i] = errB;
        // printf("%5d  %5d  %.5e  %.5e\n", m, nb, 1.0/(double)nb, errB);
    }

    // Store average into the result
    result.average = B_ave;
    // Find out the biggest error inside errors
    for(int i = 1; i <= bmax2; i = i+1){
        if(result.error < errors[i]){
            result.error = errors[i];
            result.m = pow(2, i);
        }
    }

    return result;
}