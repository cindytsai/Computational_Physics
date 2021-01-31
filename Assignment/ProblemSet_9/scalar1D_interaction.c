#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>

void exact(int N,double L,double mass,double x,double y,double *S,double *P);
void metro(int mtimes, double phimax, double *phi, double *newphi, int N, double a, double cons, double *accept, double lambda, gsl_rng *rng);
void hmc(int hmcs, double dt, int N, double a, double cons, double *phi, double *newphi, double *ph, double *accept, double lambda, gsl_rng *rng);
void ph_u(double dt, int N, double a, double cons, double *ph, double *newphi, double lambda);
void phi_u(double dt, int N, double *ph, double *newphi);

int main(){
	// Settiings
	int maxSite = 1000;
	double *phi;            // Scalar fields
	double *newphi;			// New scalar fields
	double *ph;             // Conjugate momenta.
	double a;				// lattice spacing
	double cons;			// cons = ( 2 + ma^2 ) / 2
	double accept;			// Acceptance in a single sweep
	int N;					// # of sites
	int nt;					// # of sweeps (traj) for thermalization
	int nm;					// # of sweeps (traj) for measurement
	int im;					// # of sweeps (traj) between successive measurements
	int nd;					// # of sweeps (traj) between displaying results
	int sweeps;				// Total # of sweeps (trajs) at each mass value
	int ntp;				// Total # of mass values
	int alg;				// Algorithm (0:HMC, 1:Metropolis) 
	int fw, bw;				// Forward and backward site index
	//int i, j, it;           // Variables for do loops
	int ix, iy;				// Indices for the end-points of the propagator
	int id;					// ix-iy
	int hmcs;				// # of leapfrog steps per trajectory in H.M.C
	int mtimes;				// # of updates per site in modified Metropolos

	double dt;				// d(tau) in H.M.C
	double x, y, x_end;		// end-points of the propagator.
	double L;				// length of the lattice
	double po, pn;			// Old and new field values at the site
	double so, sn;			// Old and new action of the site
	double ps;				// Sum of the neighboring field
	double phimax;			// Maximum length of phi
	double mass;			// Mass of the scalar field
	double ma;				// mass * a
	double ds;				// The change of action
	double ts;				// Accumulator for subtracted propagator
	double tp;				// Accumulator for propagator
	double tac;				// Accumulator for acceptance
	double count;			// Counter for # of measuments
	double icc;				// Counter for # of trial measurements
	double prop;			// Propagator
	double props;			// Substracted propagator
	double S, P;			// <S> and <P>
	double S_ex, P_ex;		// Exact solution of S and P
	double AC;				// Acceptance percentage
	
	double lambda;          // interaction

	gsl_rng *rng;           // Random number generator
	FILE *output;           // For output the result
	char output_name[50];   // File name
	char temp_str[10];
	
	// Initialize
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 1000);

	L = 1000.0;
	// printf("Enter the length of the lattice : \n");
	// scanf("%lf", &L);
	// printf("%lf\n", L);
	N = 100;
	// printf("Enter the # of site of the 1D lattice ( < %d ) : \n", maxSite);
	// scanf("%d", &N);
	// printf("%d\n", N);
	a = L / (double) N;

	phimax = 100.0;
	// printf("Enter the maximum length of phi : \n");
	// scanf("%lf", &phimax);
	// printf("%lf\n", phimax);
	nt = 1000;
	// printf("Enter the # of sweeps for thermalization : \n");
	// scanf("%d", &nt);
	// printf("%d\n", nt);
	nm = 1024;
	// printf("Enter # of measurements : \n");
	// scanf("%d", &nm);
	// printf("%d\n", nm);
	im = 10;
	// printf("Enter the # of sweeps between successive measurements : \n");
	// scanf("%d", &im);
	// printf("%d\n", im);
	sweeps = nt + nm * im;
	nd = 1000;
	// printf("Enter the display interval : \n");
	// scanf("%d", &nd);
	// printf("%d\n", nd);
	mass = 0.01;
	// printf("Enter mass of scalar field : \n");
	// scanf("%lf", &mass);
	// printf("%lf\n", mass);
	y = 0.0;
	// printf("Enter Initial point y :\n");
	// scanf("%lf", &y);
	// printf("%lf\n", y);
	iy = (int) (y / a + 0.5);

	printf("Enter the value of lambda : \n");
	scanf("%lf", &lambda);
	printf("%lf\n", lambda);
	printf("Choose the algorithm : \n");
	printf("0:HMC\n");
	printf("1:Metropolis\n");
	scanf("%d", &alg);
	printf("%d\n", alg);

	printf("Enter range of the final point x we want to loop through to :\n");
	scanf("%lf", &x_end);
	printf("%lf\n", x_end);

	x = 10;
	while(x <= x_end){
		// Print to file
		strcpy(output_name, "");
		strcat(output_name, "lambda");
		sprintf(temp_str, "%d", (int) (lambda * 1000000));
		strcat(output_name, temp_str);
		strcat(output_name, "x");
		sprintf(temp_str, "%d", (int) x);
		strcat(output_name, temp_str);
		strcat(output_name, ".txt");
		output = fopen(output_name, "w");

		ix = (int) (x / a + 0.5);
	
		// Initialize the grid
		phi = (double*) malloc(maxSite * sizeof(double));
		newphi = (double*) malloc(maxSite * sizeof(double));
		ph = (double*) malloc(maxSite * sizeof(double));
		for(int i = 0; i < N; i = i+1){
			// cold start
			phi[i] = 0.0;
		}

		// Choose the algorithm
		if(alg == 0){
			hmcs = 20;
			dt = 0.05;
		}
		else{
			mtimes = 6;
		}

		// Exact Solution
		exact(N, L, mass, x, y, &S_ex, &P_ex);

		// Simulation
		printf("sweeps\tmass\t\tx-y\t\t<S>\t\t<P>\t\tS_exact\t\tP_exact\t\tAC(%%)\n");
		printf("-----------------------------------------------------\n");

		ts = 0.0;
		tp = 0.0;
		tac = 0.0;
		count = 0.0;
		icc = 0.0;

		ma = mass * a;
		cons = 0.5 * pow(ma, 2) + 1.0;

		for(int i = 0; i <= sweeps; i = i+1){
			// Algorithm
			if(alg == 0){
				hmc(hmcs, dt, N, a, cons, phi, newphi, ph, &accept, lambda, rng);
			}
			else{
				metro(mtimes, phimax, phi, newphi, N, a, cons, &accept, lambda, rng);
			}

			if((i > nt) && (i % im == 0)){
				props = 0.0;
				prop = 0.0;
				for(int j = 0; j < N; j = j+1){
					po = phi[j];
					id = (j + ix - iy) % N;
					prop = prop + po * phi[id];
					props = props + po * phi[id] - pow(po, 2);
				}
				tp = tp + prop;
				ts = ts + props;
				count = count + 1;
				tac = tac + accept;
				icc = icc + 1.0;
				fprintf(output, "%e %e\n", tp / (count * N), ts / (count * N));
			}
			if((i > nt) && (i % nd == 0)){
				P = tp / (count * N);
				S = ts / (count * N);
				AC = tac / icc;
				printf("%d\t%e\t%e\t%e\t%e\t%e\t%e\t%.2f\n", i, mass, (x - y), S, P, S_ex, P_ex, AC);
			}
		}
		fclose(output);

		// Free the array
		free(phi);
		free(newphi);
		free(ph);

		// Prepare for next loop
		x = x + 10;
	}

	
	
	return 0;
}

void exact(int N,double L,double mass,double x,double y,double *S,double *P){
	int id;
	double a, R, ma, dx;

	if(mass == 0.0){
		dx = fabs(x - y);
		*S = -0.5 * dx * (L - dx) / L;
		*P = 0.0;
	}
	else{
		a = L / (double) N;
		ma = mass * a;
		id = (int) (fabs(x - y) / a + 0.5);
		R = 1.0 + 0.5 * pow(ma, 2) - ma * sqrt(1.0 + 0.25 * pow(ma, 2));
		*P = (pow(R, N - id) + pow(R, id)) / ((1.0 - pow(R, N)) * 2.0 * mass * sqrt(1.0 + 0.25 * pow(ma, 2)));
		*S = *P - (pow(R, N) + 1) / ((1.0 - pow(R, N)) * 2.0 * mass * sqrt(1.0 + 0.25 * pow(ma, 2)));
	}
}

void metro(int mtimes, double phimax, double *phi, double *newphi, int N, double a, double cons, double *accept, double lambda, gsl_rng *rng){
	double po, pn, ps, so, sn, ds;
	int fw, bw;
	*accept = 0.0;
	for(int j = 0; j < N; j = j+1){
		fw = (j + 1) % N;
		bw = (j - 1 + N) % N;
		ps = phi[fw] + phi[bw];
		for(int k = 1; k <= mtimes; k = k+1){
			po = phi[j];
			so = -po * ps + cons * pow(po, 2) + a * lambda / 24 * pow(po, 4);
			pn = (gsl_rng_uniform(rng) - 0.5) * phimax;
			sn = -pn * ps + cons * pow(pn, 2) + a * lambda / 24 * pow(pn, 4);
			ds = (sn - so) / a;
			if((ds <= 0.0) || (gsl_rng_uniform(rng) < exp(-ds))){
				phi[j] = pn;
				*accept = *accept + 1.0;
			}
		}
	}
	*accept = *accept / (N * mtimes);
}

void hmc(int hmcs, double dt, int N, double a, double cons, double *phi, double *newphi, double *ph, double *accept, double lambda, gsl_rng *rng){
	double po, ps, oldh, newh, dh;
	int fw;
	*accept = 0.0;

	for(int j = 0; j < N; j = j+1){
		ph[j] = gsl_ran_gaussian(rng, 1);
		newphi[j] = phi[j];
	}
	
	oldh=0.0;
	for(int j = 0; j < N; j = j+1){
		po = phi[j];
		fw = (j + 1) % N;
		ps = phi[fw];
		oldh = oldh + (-po * ps + cons * pow(po, 2) + a * lambda / 24 * pow(po, 4)) / a + 0.5 * pow(ph[j], 2);
	}
	ph_u(dt/2, N, a, cons, ph, newphi, lambda);
	for(int i = 1; i < hmcs; i = i+1){
		phi_u(dt, N, ph, newphi);
		ph_u(dt, N, a, cons, ph, newphi, lambda);
	}
	phi_u(dt, N, ph, newphi);
	ph_u(dt/2, N, a, cons, ph, newphi, lambda);

	newh = 0.0;
	for(int j = 0; j < N; j = j+1){
		po = newphi[j];
		fw = (j + 1) % N;
		ps = newphi[fw];
		newh = newh + (-po * ps + cons * pow(po, 2) + a * lambda / 24 * pow(po, 4)) / a + 0.5 * pow(ph[j], 2);
	}

	dh = newh - oldh;
	if(newh == oldh){
		printf("new H = old H");
	}
	else if((dh <= 0.0) || (gsl_rng_uniform(rng) < exp(-dh))){
		for(int j = 0; j < N; j = j+1){
			phi[j] = newphi[j];
		}
		*accept = 1.0;
	}
}

void ph_u(double dt, int N, double a, double cons, double *ph, double *newphi, double lambda){
	int fw, bw;
	for(int j = 0; j < N; j = j+1){
		fw = (j + 1) % N;
		bw = (j - 1 + N) % N;
		ph[j] = ph[j] + dt * (newphi[fw] + newphi[bw] - cons * 2.0 * newphi[j] - a * lambda / 6 * pow(newphi[j], 3)) / a;
		ph[j]=ph[j]+dt*(newphi[fw]+newphi[bw]-cons*2.0*newphi[j]-a*lambda/6*pow(newphi[j],3))/a;
	}
}

void phi_u(double dt, int N, double *ph, double *newphi){
	for(int j = 0; j < N; j = j+1){
		newphi[j] = newphi[j] + dt * ph[j];
	}
}
