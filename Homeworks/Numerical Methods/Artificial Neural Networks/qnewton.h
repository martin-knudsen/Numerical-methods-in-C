#ifndef HAVE_qnewton
#define HAVE_qnewton

void numeric_gradient
(double beta(gsl_vector*), gsl_vector*x, gsl_vector*grad);

int qnewton(double beta(gsl_vector*), gsl_vector*x, double acc);

#endif