#ifndef HAVE_Newton
#define HAVE_Newton
void jacobian_num(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
				  gsl_vector* x, gsl_matrix* J, double dx);

void newton(void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
	gsl_vector* xstart,
	double dx,
	double epsilon, bool ANALYTIC
);

void newton_quadratic_backtracking(
	void f(gsl_vector* x, gsl_vector* fx, gsl_matrix* J),
	gsl_vector* x,
	double dx,
	double epsilon, bool ANALYTIC
);

#endif
