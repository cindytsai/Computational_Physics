#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Calculation Used
void neighborSiteIndex(int L, int N, int *site_ip, int *site_im, int *site_jp, int *site_jm, int *site_kp, int *site_km);
double inner_product_double(int N, double *vector1, double *vector2);
float inner_product_single(int N, float *vector1, float *vector2);
void vector_ADD_double(int N, double *output, double *vector1, double *vector2, double ratio);
void vector_ADD_single(int N, float *output, float *vector1, float *vector2, float ratio);
void vector_ADD_double_single(int N, double *output, double *vector_double, float *vector_single);
void vector_laplace3D_double(int N, double mass, double *output, double *vector);
void vector_laplace3D_single(int N, float mass, float *output, float *vector);
void CG_double(int L, int N, double *b, double *x, double mass, double target_error);
void CG_single(int L, int N, float *b, float *x, float mass, float target_error);

// Algorithm
void JacobiMethod(int L, int N, double *b, double *U0, double target_error, char *output_name);
void GaussSeidelMethod(int L, int N, double *b, double *U0, double target_error, char *output_name);
void SOR(int L, int N, double *b, double *U0, double target_error, char *output_name);
void CG_Algorithm_double(int L, int N, double *b, double *x, double mass, double target_error, char *output_name);
void CG_restart(int L, int N, double *b, double *x, double mass, double target_error, double eps, int mixed_precision, char *output_name);

// Settings
// site_{ip, im, jp, jm, kp, km} -> Index of the neighboring site
int *site_ip, *site_im, *site_jp, *site_jm, *site_kp, *site_km;

int main(){
	// Settings
	// L                 -> 3D lattice box length
	// N                 -> Total number of element in array
	// site              -> Record index of the site in array
	// algorithm, mixed  -> Algorithm we want to use
	// b                 -> Source vector
	// U0                -> Matrix phi, our initial solution guess
	// target_error      -> Result must be smaller than this error
	// mass              -> For CG Algorithm, photon mass
	// output_name       -> Output file name
	int L, N, site;
	int algorithm, mixed_precision;
	double *b, *U0;
	double target_error, eps;
	double mass;
	char output_name[50];

	printf("Enter the size L of the cube L x L x L: \n");
	scanf("%d",&L);
	printf("%d\n", L);
	N = pow(L, 3);
	printf("Choose the algorithm: \n");
  	printf("0: Jacobi method\n");
  	printf("1: Gauss Seidel method\n");
  	printf("2: Successive Over-Relaxation method\n");
  	printf("3: CG Algorithm\n");
  	printf("4: CG Algorithm with restart\n"); 
  	scanf("%d", &algorithm);
  	printf("%d\n", algorithm);

  	b = (double*) malloc(N * sizeof(double));
  	U0 = (double*) malloc(N * sizeof(double));


	// Initialize the Poisson equation we want to solve
	for(int i = 0; i < L; i = i+1){
		for(int j = 0; j < L; j = j+1){
			for(int k = 0; k < L; k = k+1){
				site = i + j * L + k * pow(L, 2);
				b[site] = 0.0;
				U0[site] = 0.0;
			}
		}
	}
	b[0] = 1.0;

	// Initialize neighboring site index
	site_ip = (int*) malloc(N * sizeof(int));
	site_im = (int*) malloc(N * sizeof(int));
	site_jp = (int*) malloc(N * sizeof(int));
	site_jm = (int*) malloc(N * sizeof(int));
	site_kp = (int*) malloc(N * sizeof(int));
	site_km = (int*) malloc(N * sizeof(int));
	neighborSiteIndex(L, N, site_ip, site_im, site_jp, site_jm, site_kp, site_km);

	// Call the algorithm
	if(algorithm >= 0 && algorithm <= 2){
		printf("Input target error : \n");
		scanf("%lf", &target_error);
		printf("%.5e\n", target_error);
		printf("Input output file name: \n");
		scanf("%s", output_name);
		printf("%s\n", output_name);
		if(algorithm == 0){
			JacobiMethod(L, N, b, U0, target_error, output_name);
		}
		if(algorithm == 1){
			GaussSeidelMethod(L, N, b, U0, target_error, output_name);
		}
		if(algorithm == 2){
			SOR(L, N, b, U0, target_error, output_name);
		}
	}
	if(algorithm == 3){
		printf("Input photon mass : \n");
		scanf("%lf", &mass);
		printf("%.5e\n", mass);
		printf("Input target error : \n");
		scanf("%lf", &target_error);
		printf("%.5e\n", target_error);
		printf("Input output file name: \n");
		scanf("%s", output_name);
		printf("%s\n", output_name);
		CG_Algorithm_double(L, N, b, U0, mass, target_error, output_name);
	}
	if(algorithm == 4){
		printf("0: Double Precision (double) with Restart\n");
		printf("1: Mixed Precision (single) with Restart\n");
		scanf("%d", &mixed_precision);
		printf("%d\n", mixed_precision);
		printf("Input photon mass : \n");
		scanf("%lf", &mass);
		printf("%.5e\n", mass);
		printf("Input target error : \n");
		scanf("%lf", &target_error);
		printf("%.5e\n", target_error);
		printf("Input inner loop error : \n");
		scanf("%lf", &eps);
		printf("%.5e\n", eps);
		printf("Input output file name: \n");
		scanf("%s", output_name);
		printf("%s\n", output_name);
		CG_restart(L, N, b, U0, mass, target_error, eps, mixed_precision, output_name);
	}
	else{
		printf("Out of range\n");
	}


	return 0;
}

void neighborSiteIndex(int L, int N, int *site_ip, int *site_im, int *site_jp, int *site_jm, int *site_kp, int *site_km){
	int ip, im, jp, jm, kp, km;
	int site;

	for(int i = 0; i < L; i = i+1){
		
		im = (i - 1 + L) % L;
		ip = (i + 1) % L;
		
		for(int j = 0; j < L; j = j+1){

			jm = (j - 1 + L) % L;
			jp = (j + 1) % L;

			for(int k = 0; k < L; k = k+1){

				km = (k - 1 + L) % L;
				kp = (k + 1) % L;

				site = i + j * L + k * pow(L, 2);
				site_ip[site] = ip + j * L + k * pow(L, 2);	//Right
				site_im[site] = im + j * L + k * pow(L, 2);	//Left
				site_jp[site] = i + jp * L + k * pow(L, 2);	//Forward
				site_jm[site] = i + jm * L + k * pow(L, 2);	//Backward
				site_kp[site] = i + j * L + kp * pow(L, 2);	//Top
				site_km[site] = i + j * L + km * pow(L, 2);	//Bottom
			}
		}
	}
}

/*
 Matrix Calculation
 */
double inner_product_double(int N, double *vector1, double *vector2){
	double sum = 0.0;
	for(int i = 0; i < N; i = i+1){
		sum = sum + vector1[i] * vector2[i];
	}
	return sum;
}

float inner_product_single(int N, float *vector1, float *vector2){
	double sum = 0.0;
	for(int i = 0; i < N; i = i+1){
		sum = sum + vector1[i] * vector2[i];
	}
	return sum;
}

void vector_ADD_double(int N, double *output, double *vector1, double *vector2, double ratio){
	for(int i = 0; i < N; i = i+1){
		output[i] = vector1[i] + vector2[i] * ratio;
	}
}

void vector_ADD_single(int N, float *output, float *vector1, float *vector2, float ratio){
	for(int i = 0; i < N; i = i+1){
		output[i] = vector1[i] + vector2[i] * ratio;
	}
}

void vector_ADD_double_single(int N, double *output, double *vector_double, float *vector_single){
	for(int i = 0; i < N; i = i+1){
		output[i] = vector_double[i] + (double) vector_single[i];
	}
}

void vector_laplace3D_double(int N, double mass, double *output, double *vector){
	for(int i = 0; i < N; i = i+1){
		output[i] = (pow(mass, 2) + 6.0) * vector[i] - vector[site_ip[i]] - vector[site_im[i]] - vector[site_jp[i]] - vector[site_jm[i]] - vector[site_kp[i]] - vector[site_km[i]];
	}
}

void vector_laplace3D_single(int N, float mass, float *output, float *vector){
	for(int i = 0; i < N; i = i+1){
		output[i] = (pow(mass, 2) + 6.0) * vector[i] - vector[site_ip[i]] - vector[site_im[i]] - vector[site_jp[i]] - vector[site_jm[i]] - vector[site_kp[i]] - vector[site_km[i]];
	}
}


void CG_double(int L, int N, double *b, double *x, double mass, double target_error){
	// Solve A.x = b, and store solutions in array x
	// 
	// Settings
	// r[N]            -> Residual Vector
	// p[N]            -> Correction Vector
	// rr, bb, tt      -> (r_k, r_k), (b, b), (r_k+1, r_k+1)
	// alpha, beta     -> alpha_k, beta_k+1
	// iter            -> Record the iteration
	// max_iter        -> Maximum iteration
	// cr              -> |r_k+1| / |b|
	// cr1             -> error_k+1 - error_k

	double *r, *p, *Ap;
	double rr, bb, tt;
	double alpha, beta;
	int iter;
	int max_iter = 6000;
	double cr, cr1;

	// Initialize
	// x_0 = 0.0
	// r_0 = b - A.x_0
	// p_0 = r_0
	r = (double*) malloc(N * sizeof(double));
	p = (double*) malloc(N * sizeof(double));
	Ap = (double*) malloc(N * sizeof(double));
	for(int i = 0; i < N; i = i+1){
		x[i] = 0.0;
		r[i] = b[i];
		p[i] = r[i];
	}
	
	iter = 0;
	rr = inner_product_double(N, r, r);
	bb = inner_product_double(N, b, b);
	cr  = sqrt(rr / bb);
	cr1 = -1.0;

	while(cr > target_error && iter <= max_iter && cr1 < 0){
		// Calculate r_k+1 = r_k - alpha_k.(A.p)
		vector_laplace3D_double(N, mass, Ap, p);
		alpha = rr / inner_product_double(N, p, Ap);
		vector_ADD_double(N, r, r, Ap, -alpha);

		// Calculate x_k+1 = x_k + alpha_k.p_k
		vector_ADD_double(N, x, x, p, alpha);

		// Calculate p_k+1 = r_k+1 + beta_k+1.p_k, beta_k+1 = (r_k+1, r_k+1) / (r_k, r_k)
		tt = inner_product_double(N, r, r);
		beta = tt / rr;
		vector_ADD_double(N, p, r, p, beta);

		rr = tt;
		cr1 = cr;
		cr = sqrt(rr / bb);
		cr1 = cr - cr1;

		iter = iter + 1;
	}

	free(r);
	free(p);
	free(Ap);
}

void CG_single(int L, int N, float *b, float *x, float mass, float target_error){
	// Every thing is same as above function, but with data type "float"
	float *r, *p, *Ap;
	float rr, bb, tt;
	float alpha, beta;
	int iter;
	int max_iter = 6000;
	float cr, cr1;

	// Initialize
	// x_0 = 0.0
	// r_0 = b - A.x_0
	// p_0 = r_0
	r = (float*) malloc(N * sizeof(float));
	p = (float*) malloc(N * sizeof(float));
	Ap = (float*) malloc(N * sizeof(float));
	for(int i = 0; i < N; i = i+1){
		x[i] = 0.0;
		r[i] = b[i];
		p[i] = r[i];
	}
	
	iter = 0;
	rr = inner_product_single(N, r, r);
	bb = inner_product_single(N, b, b);
	cr  = sqrt(rr / bb);
	cr1 = -1.0;

	while(cr > target_error && iter <= max_iter && cr1 < 0){
		// Calculate r_k+1 = r_k - alpha_k.(A.p)
		vector_laplace3D_single(N, mass, Ap, p);
		alpha = rr / inner_product_single(N, p, Ap);
		vector_ADD_single(N, r, r, Ap, -alpha);

		// Calculate x_k+1 = x_k + alpha_k.p_k
		vector_ADD_single(N, x, x, p, alpha);

		// Calculate p_k+1 = r_k+1 + beta_k+1.p_k, beta_k+1 = (r_k+1, r_k+1) / (r_k, r_k)
		tt = inner_product_single(N, r, r);
		beta = tt / rr;
		vector_ADD_single(N, p, r, p, beta);

		rr = tt;
		cr1 = cr;
		cr = sqrt(rr / bb);
		cr1 = cr - cr1;

		iter = iter + 1;
	}

	free(r);
	free(p);
	free(Ap);
}

/*
Algorithm
 */
void JacobiMethod(int L, int N, double *b, double *U0, double target_error, char *output_name){
	// Settings
	// U                             -> The calculate result
	// error                         -> calculate error
	// iter, max_iter                -> Record iteration, and iteration shoud not get over max_iter
	// t_start, t_end                -> Record start of and the end of algorithm
	// time_used                     -> Calculate iteration time used 
	double *U;
	double error;
	int iter, max_iter = 6000;
	clock_t t_start, t_end;
	double time_used;
	FILE *output;

	// Initialize U
	U = (double*) malloc(N * sizeof(double));
	// for(int i = 0; i < L; i = i+1){
	// 	for(int j = 0; j < L; j = j+1){
	// 		for(int k = 0; k < L; k = k+1){
	// 			site = i + j * L + k * pow(L, 2);
	// 			U[site] = 0.0;
	// 		}
	// 	}
	// }

	// Jacobi Method, and start the clock
	t_start = clock();

	for(iter = 1; iter <= max_iter; iter = iter+1){
		// Calculate U and move result from U -> U0, for an iteration
		for(int i = 0; i < N; i = i+1){
			U[i] = (1.0/6.0) * (b[i] + U0[site_ip[i]] + U0[site_im[i]] + U0[site_jp[i]] + U0[site_jm[i]] + U0[site_kp[i]] + U0[site_km[i]]);
		}
		memcpy(U0, U, sizeof(double) * N);

		// Check the error
		error = 0.0;
		for(int i = 0; i < N; i = i+1){
			error = error + pow((U0[site_ip[i]] + U0[site_im[i]] + U0[site_jp[i]] + U0[site_jm[i]] + U0[site_kp[i]] + U0[site_km[i]]) - 6.0 * U0[i] + b[i], 2);
		}
		error = sqrt(error);
		// TODO
		printf("%d    %.5e\n", iter, error);

		if(error < target_error){
			break;
		}
		if(iter == max_iter && error > target_error){
			printf("Iteration can't reach the target error ! \n");
		}
	}

	t_end = clock();
	time_used = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;

	// Print out the result
	printf("Iteration Used: %d,  Error: %.5e\n", iter, error);
	printf("Time Used: %.10lf sec", time_used);
	output = fopen(output_name, "a");
	for(int i = 0; i < N; i = i+1){
		fprintf(output, "%.10e\n", U0[i] - U0[1]);
	}
}

void GaussSeidelMethod(int L, int N, double *b, double *U0, double target_error, char *output_name){
	// Settings
	// U                             -> The calculate result
	// error                         -> calculate error
	// iter, max_iter                -> Record iteration, and iteration shoud not get over max_iter
	// t_start, t_end                -> Record start of and the end of algorithm
	// time_used                     -> Calculate iteration time used 
	double *U;
	double error;
	int iter, max_iter = 6000;
	clock_t t_start, t_end;
	double time_used;
	FILE *output;

	// Initialize U
	U = (double*) malloc(N * sizeof(double));
	memcpy(U, U0, sizeof(double) * N);

	// Gauss-Seidel Method, and start the clock
	t_start = clock();

	for(iter = 1; iter <= max_iter; iter = iter+1){
		// Calculate U
		for(int i = 0; i < N; i = i+1){
			U[i] = (1.0/6.0) * (b[i] + U0[site_ip[i]] + U[site_im[i]] + U0[site_jp[i]] + U[site_jm[i]] + U0[site_kp[i]] + U[site_km[i]]);
		}
		memcpy(U0, U, sizeof(double) * N);

		// Check the error
		error = 0.0;
		for(int i = 0; i < N; i = i+1){
			error = error + pow((U0[site_ip[i]] + U0[site_im[i]] + U0[site_jp[i]] + U0[site_jm[i]] + U0[site_kp[i]] + U0[site_km[i]]) - 6.0 * U0[i] + b[i], 2);
		}
		error = sqrt(error);
		// TODO
		printf("%d    %.5e\n", iter, error);

		if(error < target_error){
			break;
		}
		if(iter == max_iter && error > target_error){
			printf("Iteration can't reach the target error ! \n");
		}
	}

	t_end = clock();
	time_used = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;

	// Print out the result
	printf("Iteration Used: %d,  Error: %.5e\n", iter, error);
	printf("Time Used: %.10lf sec", time_used);
	output = fopen(output_name, "a");
	for(int i = 0; i < N; i = i+1){
		fprintf(output, "%.10e\n", U0[i] - U0[1]);
	}
}

void SOR(int L, int N, double *b, double *U0, double target_error, char *output_name){
	// Settings
	// w         					 -> Relaxation factor
	// U                             -> The calculate result
	// error                         -> calculate error
	// iter, max_iter                -> Record iteration, and iteration shoud not get over max_iter
	// t_start, t_end                -> Record start of and the end of algorithm
	// time_used                     -> Calculate iteration time used 
	double w = 1.5;
	double *U;
	double error;
	int iter, max_iter = 6000;
	clock_t t_start, t_end;
	double time_used;
	FILE *output;

	// Initialize U
	U = (double*) malloc(N * sizeof(double));
	memcpy(U, U0, sizeof(double) * N);

	// Successive Over-Relation Method
	t_start = clock();

	for(iter = 1; iter <= max_iter; iter = iter+1){
		// Calculate U
		for(int i = 0; i < N; i = i+1){
			U[i] = (w/6.0) * (b[i] + U0[site_ip[i]]	+ U[site_im[i]] + U0[site_jp[i]] + U[site_jm[i]] + U0[site_kp[i]] + U[site_km[i]]) + (1.0 - w) * U0[i];
		}
		memcpy(U0, U, sizeof(double) * N);

		// Check the error
		error = 0.0;
		for(int i = 0; i < N; i = i+1){
			error = error + pow((U0[site_ip[i]] + U0[site_im[i]] + U0[site_jp[i]] + U0[site_jm[i]] + U0[site_kp[i]] + U0[site_km[i]]) - 6.0 * U0[i] + b[i], 2);
		}
		error = sqrt(error);
		// TODO
		printf("%d    %.5e\n", iter, error);

		if(error < target_error){
			break;
		}
		if(iter == max_iter && error > target_error){
			printf("Iteration can't reach the target error ! \n");
		}
	}

	t_end = clock();
	time_used = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;

	// Print out the result
	printf("Iteration Used: %d,  Error: %.5e\n", iter, error);
	printf("Time Used: %.10lf sec", time_used);
	output = fopen(output_name, "a");
	for(int i = 0; i < N; i = i+1){
		fprintf(output, "%.10e\n", U0[i] - U0[1]);
	}
}

void CG_Algorithm_double(int L, int N, double *b, double *x, double mass, double target_error, char *output_name){
	// Settings
	// r[N]            -> Residual Vector
	// p[N]            -> Correction Vector
	// rr, bb, tt      -> (r_k, r_k), (b, b), (r_k+1, r_k+1)
	// alpha, beta     -> alpha_k, beta_k+1
	// iter            -> Record the iteration
	// max_iter        -> Maximum iteration
	// cr              -> |r_k+1| / |b|
	// cr1             -> error_k+1 - error_k

	double *r, *p, *Ap;
	double rr, bb, tt;
	double alpha, beta;
	int iter;
	int max_iter = 6000;
	double cr, cr1;

	FILE *output;

	clock_t t_start, t_end;
	double time_used;

	t_start = clock();

	// Initialize
	// x_0 = 0.0
	// r_0 = b - A.x_0
	// p_0 = r_0
	r = (double*) malloc(N * sizeof(double));
	p = (double*) malloc(N * sizeof(double));
	Ap = (double*) malloc(N * sizeof(double));
	for(int i = 0; i < N; i = i+1){
		x[i] = 0.0;
		r[i] = b[i];
		p[i] = r[i];
	}
	
	iter = 0;
	rr = inner_product_double(N, r, r);
	bb = inner_product_double(N, b, b);
	cr  = sqrt(rr / bb);
	cr1 = -1.0;

	while(cr > target_error && iter <= max_iter){
		// Calculate r_k+1 = r_k - alpha_k.(A.p)
		vector_laplace3D_double(N, mass, Ap, p);
		alpha = rr / inner_product_double(N, p, Ap);
		vector_ADD_double(N, r, r, Ap, -alpha);

		// Calculate x_k+1 = x_k + alpha_k.p_k
		vector_ADD_double(N, x, x, p, alpha);

		// Calculate p_k+1 = r_k+1 + beta_k+1.p_k, beta_k+1 = (r_k+1, r_k+1) / (r_k, r_k)
		tt = inner_product_double(N, r, r);
		beta = tt / rr;
		vector_ADD_double(N, p, r, p, beta);

		rr = tt;
		cr1 = cr;
		cr = sqrt(rr / bb);
		cr1 = cr - cr1;

		iter = iter + 1;

		// TODO
		printf("%d    %.5e    %.5e\n", iter, cr, cr1);
	}

	t_end = clock();
	time_used = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;

	// Print out the result
	printf("Iteration Used: %d,  Error: %.5e\n", iter, cr);
	printf("Time Used: %.10lf sec", time_used);
	output = fopen(output_name, "a");
	for(int i = 0; i < N; i = i+1){
		fprintf(output, "%.10e\n", x[i] - x[1]);
	}

	free(r);
	free(p);
	free(Ap);
}

void CG_restart(int L, int N, double *b, double *x, double mass, double target_error, double eps, int mixed_precision, char *output_name){
	// Settings
	// Outer Loop
	// r[N]            -> Residual Vector
	// rr, bb, tt      -> (r_k, r_k), (b, b), (r_k+1, r_k+1)
	// iter            -> Record the iteration
	// max_iter        -> Maximum iteration
	// cr              -> |r_k+1| / |b|
	// cr1             -> error_k+1 - error_k
	// 
	// Inner Loop
	// t_single, t_double  -> Store Result from inner loop

	double *r, *Ax;
	double rr, bb;
	int iter;
	int max_iter = 6000;
	double cr, cr1;

	double *t_double;
	float  *t_single, *r_single;

	FILE *output;

	clock_t t_start, t_end;
	double time_used;

	t_start = clock();

	// Initialize
	// r_0 = b - A.x_0
	r = (double*) malloc(N * sizeof(double));
	Ax = (double*) malloc(N * sizeof(double));
	t_double = (double*) malloc(N * sizeof(double));
	
	t_single = (float*) malloc(N * sizeof(float));
	r_single = (float*) malloc(N * sizeof(float));
	
	// r_0 = b - A.x_0
	for(int i = 0; i < N; i = i+1){
		x[i] = 0.0;
	}
	vector_laplace3D_double(N, mass, Ax, x);
	vector_ADD_double(N, r, b, Ax, -1.0);
	
	iter = 0;
	rr = inner_product_double(N, r, r);
	bb = inner_product_double(N, b, b);
	cr  = sqrt(rr / bb);
	cr1 = -1.0;

	while(cr > target_error && iter <= max_iter){
		
		// CG restart with mixed precision (single)
		if(mixed_precision == 1){
			// Assign r (double) to r_single (float)
			for(int i = 0; i < N; i = i+1){
				r_single[i] = (float) r[i];
			}
			// Inner Loop with single precision
			CG_single(L, N, r_single, t_single, (float)mass, (float)eps);

			// Calculate x_k+1 = x_k + t_single_k, note: x is double precision
			vector_ADD_double_single(N, x, x, t_single);
		}

		// CG restart with double precision (double)
		if(mixed_precision == 0){
			CG_double(L, N, r, t_double, mass, eps);
			vector_ADD_double(N, x, x, t_double, 1.0);
		}

		vector_laplace3D_double(N, mass, Ax, x);
		vector_ADD_double(N, r, b, Ax, -1.0);

		// Calculate error
		cr1 = cr;
		rr = inner_product_double(N, r, r);
		cr = sqrt(rr / bb);
		cr1 = cr - cr1;

		iter = iter + 1;

		// TODO
		printf("%d    %.5e    %.5e\n", iter, cr, cr1);
	}

	t_end = clock();
	time_used = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;

	// Print out the result
	printf("Iteration Used: %d,  Error: %.5e\n", iter, cr);
	printf("Time Used: %.10lf sec", time_used);
	output = fopen(output_name, "a");
	for(int i = 0; i < N; i = i+1){
		fprintf(output, "%.10e\n", x[i] - x[1]);
	}

	free(r);
	free(Ax);
	free(t_double);
	free(t_single);
	free(r_single);
}
