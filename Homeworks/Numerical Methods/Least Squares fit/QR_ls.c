#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

void least_squares(gsl_vector* x, gsl_vector* y, gsl_vector* dy, int m,
				 double funs(int i, double x), gsl_vector* c, gsl_matrix* COV){
	int n = x->size;

	gsl_matrix* A = gsl_matrix_alloc(n, m);
	gsl_matrix* AT = gsl_matrix_alloc(m, n);
	gsl_matrix* A_clone = gsl_matrix_alloc(n, m);
	gsl_matrix* R = gsl_matrix_alloc(m, m);
	gsl_matrix* RINV = gsl_matrix_alloc(m, m);
	gsl_matrix* ATA = gsl_matrix_alloc(m, m);
	gsl_matrix* I = gsl_matrix_alloc(m, m);
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
			gsl_matrix_set(AT,k,i,Aik);
		}
	}

	qr_gs_decomp(A,R);
	// now Q=A
    qr_gs_solve(A,R,b,c);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AT, A, 0.0, ATA);

    /*
    Here comes a cool trick, which I must admit I learned from dmitris code
    We need to find COV using (4.14). A smart way to do this would be
    if we could somehow find R^-1. We could use our function for calculating
    the inverse using gram-schmidt. However that function takes the 
    QR-decomp of whatever we are finding the the inverse of. 
    So what is it of R? Well we already know that R is upper 
    triangular, which is the only thing we need R to be. 
    So what is the orthogonal matrix Q? how about just the identity
    matrix! Now we have that R=IR, the product of an orthogonal matrix
    and an upper triangular one and voila!
     */ 

    gsl_matrix_set_identity(I);
    qr_gs_inverse(I, R, RINV);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, RINV, RINV, 0.0, COV);

	gsl_matrix_free(A);
	gsl_matrix_free(A_clone);
	gsl_matrix_free(R);
	gsl_matrix_free(ATA);
	gsl_matrix_free(AT);
	gsl_matrix_free(RINV);
	gsl_matrix_free(I);
	gsl_vector_free(b);
	gsl_vector_free(b_clone);
} 
	