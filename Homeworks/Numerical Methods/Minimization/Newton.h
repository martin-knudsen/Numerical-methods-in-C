#ifndef HAVE_Newton
#define HAVE_Newton

void newton_minimization(double f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
	gsl_vector* xstart,
	double dx,
	double epsilon
);

#endif
