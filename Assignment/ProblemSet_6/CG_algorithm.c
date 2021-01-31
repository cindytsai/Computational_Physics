#include<stdio.h>
#include<tgmath.h>
#include<math.h>
#include<complex.h>

double complex inner_product(int N, double complex *vector1, double complex *vector2);
void matrix_Dagger(int N, double complex *output, double complex *matrix);
void vector_ADD(int N, double complex *output, double complex *vector1, double complex *vector2, double complex ratio);
void vector_MUL(int N, double complex *output, double complex *matrix, double complex *vector);
void CG_Algorithm(int N, double complex *D, double complex *b);

int main(void){
	// Settings
	// N            -> Length of the square matrix D
	// D[N*N]       -> Store the Matrix Dx=b as a 1D double complex array
	// b[N]         -> Store the vector 1D double complex array
	// inputRe      -> input file , Real part
	// inputIm      -> input file , Imagine part
	// temp1, temp2 -> store input temperary
	int N = 5;
	double complex D[N*N], b[N];
	FILE *inputRe, *inputIm;
	double temp1, temp2;

	// Testing CG with other matrix
	double complex A[25], c[5];

	// Input Complex Matrix D
	inputRe = fopen("D_Re.txt", "r");
	inputIm = fopen("D_Im.txt", "r");
	for(int i = 0; i < N*N; i = i+1){
		fscanf(inputRe, "%lf", &temp1);
		fscanf(inputIm, "%lf", &temp2);
		D[i] = temp1 + temp2 * I;
	}
	fclose(inputRe);
	fclose(inputIm);

	// Input Complex Vector b
	inputRe = fopen("b_Re.txt", "r");
	inputIm = fopen("b_Im.txt", "r");
	for(int i = 0; i < N; i = i+1){
		fscanf(inputRe, "%lf", &temp1);
		fscanf(inputIm, "%lf", &temp2);
		b[i] = temp1 + temp2 * I;
	}
	fclose(inputRe);
	fclose(inputIm);

	// Print Complex Matrix D
	printf("D = \n");
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			printf("%.5e+%.5e*I ", creal(D[i*N+j]), cimag(D[i*N+j]));
		}
		printf("\n");
	}
	// Print Complex Vector b
	printf("b = \n");
	for(int i = 0; i < N; i = i+1){
		printf("%.5e+%.5e*I\n", creal(b[i]), cimag(b[i]));
	}
	
	// Call CG_Algorithm
	printf("\nSolve Dx = b\n");
	CG_Algorithm(N, D, b);

	// // Test CG_Algorithm with Matrix in Q1
	// for(int i = 0; i < 5; i = i+1){
	// 	for(int j = 1; j <= 5; j = j+1){
	// 		A[(i*5) + (j-1)] = (1.0 / (double) (i + j)) + 0.0 * I;
	// 	}
	// }
	// for(int i = 0; i < 5; i = i+1){
	// 	c[i] = 0.0 + 0.0 * I;
	// 	if(i == 0) c[i] = 1.0 + 0.0 * I;
	// }

	// printf("\nSolve Ax = c\n");
	// CG_Algorithm(N, A, c);
	
	return 0;
}

double complex inner_product(int N, double complex *vector1, double complex *vector2){
	// Compute inner product , and return the output
	// (vector1, vector2) = vector1.dagger * vector2
	double complex sum = 0.0 + 0.0 * I;
	for(int i = 0; i < N; i = i+1){
		sum = sum + conj(vector1[i]) * vector2[i];
	}
	return sum;
}

void matrix_Dagger(int N, double complex *output, double complex *matrix){
	// Add dagger to a matrix , and return the output as 1D array
	// output = matrix.dagger (transpose and conjugate)
	for(int i = 0; i < N; i = i+1){
		for(int j = 0; j < N; j = j+1){
			output[i*N+j] = conj(matrix[i+j*N]);
		}
	}
}

void vector_ADD(int N, double complex *output, double complex *vector1, double complex *vector2, double complex ratio){
	// Compute Addition of vector, and the output is vector
	// output = vector1 + ratio * vector2
	for(int i = 0; i < N; i = i+1){
		output[i] = vector1[i] + vector2[i] * ratio;
	}
}

void vector_MUL(int N, double complex *output, double complex *matrix, double complex *vector){
	// Compute Multiplication of Matrix * vector, and the output is vector
	// output = matrix * vector
	double complex sum;
	for(int i = 0; i < N; i = i+1){
		sum = 0.0 + 0.0 * I;
		for(int j = 0; j < N; j = j+1){
			sum = sum + matrix[i*N+j] * vector[j];
		}
		output[i] = sum;
	}
}

void CG_Algorithm(int N, double complex *D, double complex *b){
	// Settings
	// 
	// D[N*N]       -> Store matrix Dx = b , in 1D double complex array D
	// b[N]         -> Store matrix Dx = b , in 1D double complex array b
	// r[N]         -> Residual vector
	// s[N]         -> Gradient
	// p[N]         -> Correction vector
	// x[N]         -> Solution x
	// 
	// D_d[N*N]     -> D.dagger
	// eps          -> epsilon, accept error
	// cr           -> criteria, calculate error
	// iteration    -> count how many iteration
	// lambda,mu    -> calculate lambda, mu
	// Dp[N]        -> Store result of (D * p_k) vector
	// DDp[N]       -> Store result of (D.dagger * D * p_k) vector
	double complex r[N], s[N], p[N], x[N];
	double complex D_d[N*N];
	double eps = 1.0e-10;
	double cr;
	double complex lambda, mu, Dpk2;
	int iteration = 0;
	double complex Dp[N], DDp[N];

	// Initialize 
	// D_d = D.dagger
	// x_0 = 0, r_0 = b, s_0 = D.dagger * r_0, p_0 = s_0
	matrix_Dagger(N, D_d, D);
	for(int i = 0; i < N; i = i+1){
		x[i] = 0.0;
		r[i] = b[i];
	}
	vector_MUL(N, s, D_d, r);
	for(int i = 0; i < N; i = i+1){
		p[i] = s[i];
	}
	cr = creal(inner_product(N, r, r));

	// Start the iteration
	while(cr > eps){
		// Add one iteration
		iteration = iteration + 1;

		// Calculate lambda
		// lambda = (pk, sk) / (Dpk, Dpk)
		vector_MUL(N, Dp, D, p);
		Dpk2 = inner_product(N, Dp, Dp);
		lambda = inner_product(N, p, s);
		lambda = lambda / Dpk2;

		// Calculate x_k+1 = x_k + lambda_k * p_k
		vector_ADD(N, x, x, p, lambda);
		// Calculate r_k+1 = b - D * x_k+1 = r_k - lambda_k * D * p_k
		vector_ADD(N, r, r, Dp, (-1.0) * lambda);
		// Calculate s_k+1 = D.dagger * r_k+1;
		vector_MUL(N, s, D_d, r);

		//Calculate mu
		// mu = -(D.dagger * D * p_k , s_k+1) / (D * p_k, D * p_k)
		vector_MUL(N, DDp, D_d, Dp);
		mu = (-1.0) * inner_product(N, DDp, s) / Dpk2;

		// Calculate p_k+1 = s_k+1 + mu_k * p_k
		vector_ADD(N, p, s, p, mu);

		// Calculate criteria , cr = |r_k+1| ^ 2 = (r_k+1, r_k+1)
		cr = creal(inner_product(N, r, r));

		// Print iteration and criteria(error)
		printf("Iteration : %d ,  Error = %.5e\n", iteration, cr);
	}

	// Print the solution x
	printf("x = \n");
	for(int i = 0; i < N; i = i+1){
		printf("x%d = %.5e + (%.5e) * I\n", i+1, creal(x[i]), cimag(x[i]));
	}
	
}