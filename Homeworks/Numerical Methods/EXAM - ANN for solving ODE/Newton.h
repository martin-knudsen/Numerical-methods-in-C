#ifndef HAVE_Newton
#define HAVE_Newton

void gradient_num(void f(gsl_vector* x, gsl_vector* fx),
				  gsl_vector* x, gsl_vector* df, double dx);

void newton_minimization(double f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
	gsl_vector* xstart,
	double dx,
	double epsilon
);

void newton_minimization_Broyden(
	double f(gsl_vector* x, gsl_vector* df),
	gsl_vector* x,
	double dx,
	double epsilon
);

void newton_minimization_Broyden_num_gradient(
	double f(gsl_vector* x),
	gsl_vector* x,
	double dx,
	double epsilon
);

#endif
