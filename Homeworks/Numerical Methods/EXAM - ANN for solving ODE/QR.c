#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f"

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {
	int m = A->size2;

	gsl_vector_view a_i;
	gsl_vector_view a_j;
	
	double R_ij, R_ii;
	for(int i=0; i<m; i++){
		a_i= gsl_matrix_column(A, i);
		R_ii=gsl_blas_dnrm2(&a_i.vector);
		gsl_matrix_set(R, i, i, R_ii);
		gsl_vector_scale(&a_i.vector, 1.0/R_ii);
		
		for(int j=i+1; j<m; j++){
			a_j = gsl_matrix_column(A, j);
			gsl_blas_ddot(&a_i.vector, &a_j.vector, &R_ij);
			gsl_matrix_set(R, i, j, R_ij);
			gsl_matrix_set(R,j,i,0);
			gsl_blas_daxpy(-R_ij,&a_i.vector,&a_j.vector);
		}
	}
}

void qr_gs_solve(gsl_matrix *Q, gsl_matrix* R, \
	gsl_vector *b,gsl_vector *x){

	//gsl_blas_dgemv(CblasNoTrans, 1.0, Q, b, 0.0, x);
	gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);

	/* backsubstitution from Dmitri Fedorovs following (2.4)
		first time skipping the first for loop for starting point.
		then iterating down through yi each time first getting ci then
		subtracting the sum of Uik yik and dividing by Uii. 

	*/
	for(int i=x->size-1; i>=0; i--){
		double s = gsl_vector_get(x,i);
		for(int k=i+1; k<x->size; k++){
			s -= gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));
	}
}

void qr_solve_system(gsl_matrix* A, gsl_vector* b, gsl_vector* x){
	int n = A->size1;
	int m = A->size2;

	gsl_matrix* A_clone = gsl_matrix_alloc(n,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);

	gsl_matrix_memcpy(A_clone, A);
	qr_gs_decomp(A_clone,R);
	qr_gs_solve(A_clone,R,b,x);

	gsl_matrix_free(A_clone);
	gsl_matrix_free(R);
}

void qr_gs_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	int n = Q->size1;

	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* e = gsl_vector_alloc(n);

	for(int i=0; i<n; i++) {
		gsl_vector_set(e, i, 1.0);
		qr_gs_solve(Q,R,e,x);
		gsl_matrix_set_col(B, i, x);
		gsl_vector_set(e, i, 0.0);
	}
	gsl_vector_free(x);
	gsl_vector_free(e);
}