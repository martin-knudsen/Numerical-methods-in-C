#include <stdio.h>
#include <math.h>
#include <gsl/gsl_multimin.h>


double fit_equation(double *t, double *A, double *B, double *T) {
	double f_val = (*A)*exp(-(*t)/(*T))+(*B);

	return f_val; 
}

int master() {
	
}

int main() {



	return EXIT_SUCCESS;
}