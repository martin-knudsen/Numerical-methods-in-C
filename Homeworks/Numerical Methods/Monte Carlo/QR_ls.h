#ifndef HAVE_QR_ls
#define HAVE_QR_ls

void least_squares(gsl_vector* x, gsl_vector* y, gsl_vector* dy, int m,
				 double funs(int i, double x), gsl_vector* c, gsl_matrix* result);

#endif