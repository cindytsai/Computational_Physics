#include<stdio.h>
#include<math.h>

// dy/dx = f(x, y)
// exact_sol = y = tan(x)
double f(double x, double y);
double exact_sol(double x);
// ODE solver
double Euler(double x, double end_x, double y, double h);
double ModifiedEuler(double x, double end_x, double y, double h);
double ImprovedEuler(double x, double end_x, double y, double h);
double RungeKutta_3rd(double x, double end_x, double y, double h);

int main(void){
	// Settings
	// step          -> Record how many steps we have calculate
	// h             -> Step size
	// start, end    -> Do the evaluation from start ~ end
	// x, y[4]       -> Record x and y = {Euler, ModifiedEuler, ImprovedEuler, RungeKutta_3rd}
	int step = 1;
	double h = 0.01;
	double start = 0.0, end = 1.0;
	double x, y[4];

	// Initialized
	x = start;
	for(int i = 0; i < 4; i = i+1){
		y[i] = exact_sol(x);
	}
	
	printf("Step    x       Exact_sol    Euler      ModifiedEuler    ImprovedEuler    RungeKutta3rd\n");
	printf("______________________________________________________________________________________________\n");

	while(x < end){

		// Call the Algorithm
		y[0] = Euler(x, x+h, y[0], h);
		y[1] = ModifiedEuler(x, x+h, y[1], h);
		y[2] = ImprovedEuler(x, x+h, y[2], h);
		y[3] = RungeKutta_3rd(x, x+h, y[3], h);

		// Print the result
		printf("%d\t", step);
		printf("%.2f\t", x+h);
		printf("%.5f      ", exact_sol(x+h));
		for(int i = 0; i < 4; i = i+1){
			printf("%.5f      ", y[i]);
		}
		printf("\n");

		// Add one step
		x = x + h;
		step = step + 1;
	}

	return 0;
}

double f(double x, double y){
	return pow(y, 2) + 1;
}

double exact_sol(double x){
	return tan(x);
}

double Euler(double x, double end_x, double y, double h){
	// Settings
	// x            -> Start from x = start
	// end_x        -> End when x = end
	// y            -> Record y
	// h            -> Step size

	while(x < end_x){

		y = y + h * f(x, y);

		x = x + h;
	}
	
	return y;
}

double ModifiedEuler(double x, double end_x, double y, double h){
	// Settings
	// x            -> Start from x = start
	// end_x        -> End when x = end
	// y            -> Record y
	// h            -> Step size
	double k1, k2;
	
	while(x < end_x){

		k1 = f(x, y);
		k2 = f(x + (0.5 * h), y + (0.5 * h * k1));
		y = y + h * k2;

		x = x + h;
	}

	return y;
}

double ImprovedEuler(double x, double end_x, double y, double h){
	// Settings
	// x            -> Start from x = start
	// end_x        -> End when x = end
	// y            -> Record y
	// h            -> Step size
	double k1, k2;

	while(x < end_x){

		k1 = f(x, y);
		k2 = f(x + h, y + h * k1);
		y = y + (h / 2.0) * (k1 + k2);

		x = x + h;
	}
	
	return y;
}

double RungeKutta_3rd(double x, double end_x, double y, double h){
	// Settings
	// x            -> Start from x = start
	// end_x        -> End when x = end
	// y            -> Record y
	// h            -> Step size
	double k1, k2;

	while(x < end_x){

		k1 = f(x, y);
		k2 = f(x + (2.0/3.0) * h, y + (2.0/3.0) * h * k1);
		y = y + (h / 4.0) * (k1 + 3.0 * k2);

		x = x + h;
	}

	return y;
}