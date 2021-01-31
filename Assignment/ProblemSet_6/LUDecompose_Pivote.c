#include<stdio.h>
#include<math.h>
#include<string.h>

void LU_Decomposition(int N, double *A_1D);


int main( void ){

	// test matrix
	double A_1D[9] = {1, 1, 0, 
				   	  1, 1, 2,
				   	  1, 2, 1};
	// Target matrix
	double B_1D[25];
	for(int i = 0; i < 5; i = i+1){
		for(int j = 1; j <= 5; j = j+1){
			B_1D[(i*5) + (j-1)] = 1.0 / (double) (i + j);
		}
	}

	LU_Decomposition(3, A_1D);
	printf("==========================================\n");
	LU_Decomposition(5, B_1D);

	return 0;
}


void LU_Decomposition(int N, double *A_1D){
	// Settings
	// 
	// For LU-decomposition
	// L[N][N], U[N][N]  -> Store matrix L and U
	// A                 -> Our target
	// P[N]              -> Record permutation
	// lu_sum            -> Store sigma(L*U)
	// swap_row          -> Record row we want to swap
	// flag              -> If we encounter L[i][i] == 0 then flag = 1
	// 
	// For inverse matrix
	double L[N][N], U[N][N];
	double A[N][N];
	int P[N];
	double lu_sum;
	int swap_row, check_row, flag = 0;

	// Initialize L, U, P, and the diagonal of U is 1.0
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			L[i][j] = 0.0;
			U[i][j] = 0.0;
			if(i == j) U[i][j] = 1.0;
		}
		P[i] = i;
	}
	// Reshape matrix A_1D to A
	for(int i = 0; i < N*N; i = i+1){
		A[i/N][i%N] = A_1D[i];
	}
	
	// Print A
	printf("A = \n");
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			printf("%.5e ", A[i][j]);
		}
		printf("\n");
	}
	
	// LU-decomposition of PA
	flag = 0;
	check_row = 0;
	swap_row = check_row + 1;

	do{

		for(int k = 0; k < N; k = k+1){
			
			// Compute L
			for(int i = k; i < N; i = i+1){

				// Compute element of L[i][k]
				lu_sum = 0.0;
				for(int j = 0; j < k; j = j+1){
					lu_sum = lu_sum + L[i][j] * U[j][k];
				}
				L[i][k] = A[P[i]][k] - lu_sum;

				// Check if L[i][i] == 0.0, if so, swap row i with row below i ( swap_row > i )
				if (i == k && L[i][k] == 0.0){

					if(swap_row >= N){
						printf("Having Zero Pivote ! det(A) = 0\n");
						flag = 0;
					}
					else{
						P[i] = swap_row;
						P[swap_row] = i;
						flag = 1;
					}
					swap_row = swap_row + 1;
					
					break;
				}

				if(i == k && L[i][k] != 0.0){
					check_row = check_row + 1;
					swap_row = check_row + 1;
					flag = 0;
				}

			}

			if(flag == 1) break;

			// Compute U
			for(int j = k; j < N; j = j+1){
				lu_sum = 0.0;
				for(int i = 0; i < k; i = i+1){
					lu_sum = lu_sum + L[k][i] * U[i][j];
				}
				U[k][j] = (1.0 / L[k][k]) * (A[P[k]][j] - lu_sum);
			}
		}

	}while(flag != 0);

	// Print the result P, PA, L, U
	printf("P = \n");
	for(int i = 0; i < N; i = i+1){
		printf("%d ", P[i]);
	}
	printf("\n");

	printf("PA = \n");
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			printf("%.5e ", A[P[i]][j]);
		}
		printf("\n");
	}

	printf("L = \n");
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			printf("%.5e ", L[i][j]);
		}
		printf("\n");
	}

	printf("U = \n");
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			printf("%.5e ", U[i][j]);
		}
		printf("\n");
	}

}