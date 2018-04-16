#ifndef HAVE_QR
#define HAVE_QR

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);

void qr_gs_solve(gsl_matrix *Q, gsl_matrix* R, \
	gsl_vector *b,gsl_vector *x);

void qr_solve_system(gsl_matrix* A, gsl_vector* b, gsl_vector* x);

void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B);

#endif
