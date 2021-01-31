#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double y(double Bjk_s);
void getBinningJackknife(double *A, int bmax2, int N);

int main(){
	// Settings
	// N successive measurements
	// A[1~N] for storing data --> Simple Sampling
	// B[1~N] for storing data --> Metropolis Algorithm
    // 
	// t for computing the integrated autocorrelation time, from t : 0 ~ tmax
    // bmax2 for max of the numbers of blocks, 2^bmax2
    // 
    // Store best estimate of the secondary quantity <y> = y(<A>) --> estA, estB;
	int N = pow(2, 10);
	double *A, *B;
    int bmax2 = 10;
    double esty_A = 0.0, esty_B = 0.0;

    A = (double*) malloc(N * sizeof(double));
    B = (double*) malloc(N * sizeof(double));

	// Read from data.txt
    FILE *infile;
    infile = fopen("data.txt","r");
    for(int i = 0; i < N; i = i+1) {
      fscanf(infile,"%lf %lf", &A[i], &B[i]);
    }
    fclose(infile);

    // Find best estimate of the secondary qunatity <y> = y(<A>)
    for(int i = 0; i < N; i = i+1){
        esty_A = esty_A + A[i];
        esty_B = esty_B + B[i]; 
    }
    esty_A = y(esty_A / (double)N);
    esty_B = y(esty_B / (double)N);

    printf("Best estimate of the secondary quantity:\n");
    printf("Simple Sampling      : %.5e\n", esty_A);
    printf("Metropolis Algorithm : %.5e\n\n", esty_B);

    // Call getBinningJackknife Method
    printf("____Simple Sampling (Binning with Jackknife)____\n");
    getBinningJackknife(A, bmax2, N);
    printf("____Metropolis Algorithm (Binning with Jackknife)____\n");
    getBinningJackknife(B, bmax2, N);

    return 0;
}

double y(double Bjk_s){
    return exp(-Bjk_s);
}

void getBinningJackknife(double *A, int bmax2, int N){
    // General Settings
    // numbers of blocks m from : 2^1 ~ 2^bmax2
    // numbers of data per block --> nb
    // block average --> B[i]
    int m, nb;
    double *B, *Bjk;
    double yjk_ave, errorJK;


    printf("  m      nb       yjk_ave      errorJK\n");

    for(int i = bmax2; i >= 1; i = i-1){
        
        m = pow(2, i);
        nb = N / m;
        B = (double*) malloc(m * sizeof(double));
        Bjk = (double*) malloc(m * sizeof(double));
        yjk_ave = 0.0;
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

        // Calculate yjk_ave
        for(int j = 0; j < m; j = j+1){
            yjk_ave = yjk_ave + y(Bjk[j]);
        }
        yjk_ave = yjk_ave / (double)m;

        // Calculate errorJK
        for(int j = 0; j < m; j = j+1){
            errorJK = errorJK + pow((y(Bjk[j]) - yjk_ave), 2);
        }
        errorJK = errorJK * ((double)(m - 1) / (double)m);
        errorJK = sqrt(errorJK);

        // Print out the result
        printf("%5d  %5d  %.5e  %.5e\n", m, nb, yjk_ave, errorJK);

    }
    printf("\n");
}