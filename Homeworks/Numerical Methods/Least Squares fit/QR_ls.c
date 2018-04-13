#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void least_squares(gsl_vector* x, gsl_vector* y, gsl_vector* dy, int m,
				 double funs(int i, double x), gsl_vector* c, gsl_matrix* COV){
	int n = x->size;

	gsl_matrix* A = gsl_matrix_alloc(n, m);
	gsl_matrix* A_clone = gsl_matrix_alloc(n, m);
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_vector* b_clone = gsl_vector_alloc(n);


	double bi, xi, yi, dyi, Aik;
	for(int i=0; i<n; i++) {
		xi = gsl_vector_get(x, i);
		yi = gsl_vector_get(y, i);
		dyi = gsl_vector_get(dy, i);
		bi = yi/dyi;
		gsl_vector_set(b, i, bi);		
		for(int k=0; k<m; k++){
			Aik = funs(k, xi)/dyi;
			gsl_matrix_set(A,i,k,Aik);
			gsl_matrix_set(A_clone,i,k,Aik);
		}
	}

	qr_gs_decomp(A,R);
	// now Q=A
    qr_gs_solve(A,R,b,c);


	gsl_matrix_free(A);
	gsl_matrix_free(A_clone);
	gsl_matrix_free(R);
	gsl_vector_free(b);
	gsl_vector_free(b_clone);
} 
	