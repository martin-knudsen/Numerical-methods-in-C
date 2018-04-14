#ifndef HAVE_s_ls
#define HAVE_s_ls

void singular_decomp(gsl_matrix* A, gsl_matrix* U, 
					 gsl_matrix* S, gsl_matrix* V);

void least_squares_problem_singular(gsl_matrix* A, gsl_vector* x, 
									gsl_vector* b, gsl_matrix* COV) ;

void least_squares_singular(gsl_vector* x, gsl_vector* y, gsl_vector* dy, int m,
				 double funs(int i, double x), gsl_vector* c, gsl_matrix* COV);

#endif