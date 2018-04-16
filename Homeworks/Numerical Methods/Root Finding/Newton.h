#ifndef HAVE_Newton
#define HAVE_Newton

void newton(void f(gsl_vector* x, gsl_vector* fx),
	gsl_vector* xstart,
	double dx,
	double epsilon
);

#endif
