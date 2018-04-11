#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX;

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {
	int m = A->size2;
	int n = A->size1;

	gsl_vector* a_i = gsl_vector_alloc(n);
	gsl_vector* q_i = gsl_vector_alloc(n);
	gsl_vector* a_j = gsl_vector_alloc(n);

	int status; double dotproduct, R_ii, R_ij;
	for(int i=0; i<m; i++){
		gsl_matrix_get_col(a_i, A, i);
		status= gsl_blas_ddot(a_i, a_i, dotproduct);
		R_ii=sqrt(dotproduct);
		gsl_matrix_set(R, i, i, R_ii);
		status = int gsl_vector_memcpy(q_i, a_i);
		status = int gsl_vector_scale(q_i, 1.0/R_ii);
		status = int gsl_matrix_set_col(A, i, q_i);

		for(int j=i+1; j<m; i++){
			gsl_matrix_get_col(a_j, A, j);
			status= gsl_blas_ddot(q_i, a_j, R_ij);
			gsl_matrix_set(R, i, j, R_ij);
			status = int gsl_vector_scale(q_i, R_ij);
			status = int gsl_vector_sub(a_j, q_i);
			status = gsl_matrix_set_col(A, j, a_j);
		}

	}




	gsl_vector_free(a_i);
	gsl_vector_free(a_j);
	gsl_vector_free(q_i);
	
}
/*
void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, \
	const gsl_vector* b,gsl_vector* x){

}
*/
int main() {
	
	const int n = 10;
	const int m = 5

	gsl_matrix* A = gsl_matrix_alloc(n, m);
	gsl_matrix* R = gsl_matrix_alloc(m, m);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			gsl_matrix_set(A, i, j, RND);
		}
	}

	qr_gs_decomp(A, R);

	printf("checking that R is upper triangular\n");
	printf("by checking wether left-down is zero and rest isn't:\n");

	double R_ij; bool status=true;
	for(int i=1; i<n; i++){
		for(int j=0; j<i; j++){
			R_ij = gsl_matrix_get(R, i, j);
			if(R_ij!=0) status=false;
		}
	}

	printf(status ? "True\n": "False\n");

	printf("checking Q^TQ=1 by using gsl_BLAS\n:");

	gsl_matrix* A_T = gsl_matrix_alloc(m,n);

	status_int = gsl_matrix_transpose_memcpy(A_T, A);

	gsl_matrix* C = gsl_matrix_alloc(m, m);

	status_int = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_T, A, 0.0, C);

	


	return EXIT_SUCCESS;
}