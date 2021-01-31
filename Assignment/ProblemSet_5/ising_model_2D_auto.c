#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

struct Correlation{
    int t;
    double tauA;
};

struct Correlation getAutocorrelationTime(char *inputfile, int N, int tmax){
    // Settings
    // inputfile      -> name of the input file
    // N              -> number of measurements, read how many line of input file
    // tmax           -> computing the integrated autocorrelation time, from t : 0 ~ tmax
    // A_ave, A2_ave  -> <A> <A^2>
    // A[1~N]         -> Store data from input file
    // t, tauA        -> t, tauA
    // sat            -> saturation or not, if reach saturation sat = 1
    int num;
    double *A, *A_cor;
    double A_ave = 0.0, A2_ave = 0.0;
    double C;
    struct Correlation ans = {.t = 0, .tauA = 0.0};
    int sat;
    double temp;

    A = (double*) malloc((N + 1) * sizeof(double));
    A_cor = (double*) malloc((tmax + 1) * sizeof(double));

    // Read from file name inputfile
    FILE *infile;
    infile = fopen(inputfile,"r");
    for(int i = 1; i <= N; i = i+1) {
      fscanf(infile,"%lf %lf", &temp, &A[i]);
    }
    fclose(infile);

    // Calculate <A>, <A^2>
    for(int i = 1; i <= N; i = i+1){
        A_ave = A_ave + A[i];
        A2_ave = A2_ave + pow(A[i], 2);
    }
    A_ave = A_ave / (double)N;
    A2_ave = A2_ave / (double)N;

    // Calculate C(j)
    for(int j = 0; j <= tmax; j = j+1){
        A_cor[j] = 0.0;
        num = 0;
        for(int i = 1; i <= N-j; i = i+1){
            A_cor[j] = A_cor[j] + A[i] * A[i+j];
            num = num + 1;
        }
        A_cor[j] = (A_cor[j] / (double)num) - pow(A_ave, 2);
    }

    // Find integrated autocorrelation time
    sat = 0;
    for(int t = 1; t <= tmax ; t = t+1){
        C = A_cor[t] / A_cor[0];

        // Reach saturation
        if( C < 0.0){
            ans.t = t - 1;
            sat = 1;
            break;
        }

        ans.tauA = ans.tauA + C;
    }
    // If doesn't reach saturation
    if(sat != 1){
        printf("File Name : %s, not saturate !\n", inputfile);
    }

    return ans;
}

int main(int argc, char *argv[]){

    char filename[100];
    int nm = 10000;
    int tmax = 10000;
    struct Correlation output;

    strcpy(filename, argv[1]);
    output = getAutocorrelationTime(filename, nm, tmax);

    // printf("  t    tauA\n");
    printf("%d %.5e\n", output.t, output.tauA);

    return 0;
}