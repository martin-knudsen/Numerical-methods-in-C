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

// inspired by Dmitri Fedorovs print implementation
void printm(gsl_matrix *A){
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) {
			printf(FMT,gsl_matrix_get(A,i,j));
		}
		printf("\n");}
}

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R) {
	int m = A->size2;
	int n = A->size1;

	gsl_vector_view a_i;
	gsl_vector_view a_j;
	
	int status; double R_ij, R_ii;
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

void qr_gs_solve(const gsl_matrix* Q, const gsl_matrix* R, \
	const gsl_vector* b,gsl_vector* x){

	gsl_blas_dgemv(CblasTrans, 1.0, Q, b, 0.0, x);

	for(int i=x−>size−1; i>=0; i−−){
		double s = gsl_vector_get(x,i);
		for(int k=i+1; k<x−>size; k++){
			s −= gsl_matrix_get(R,i,k)∗gsl_vector_get(x,k);
		}
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}
}

int main() {
	
	const int n = 6;
	const int m = 4;

	gsl_matrix* A = gsl_matrix_calloc(n, m);
	gsl_matrix* A_clone = gsl_matrix_calloc(n, m);
	gsl_matrix* R = gsl_matrix_calloc(m, m);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			gsl_matrix_set(A, i, j, RND);
		}
	}

	int status_int = gsl_matrix_memcpy(A_clone, A);
	qr_gs_decomp(A, R);

	printf("This is A:\n");
	printm(A_clone);
	printf("This is Q:\n");
	printm(A);
	printf("This is R:\n");
	printm(R);
	printf("Checking that R is upper triangular\n");
	printf("By checking wether left-down is zero and rest isn't:\n");

	double R_ij; bool status=true;
	for(int i=1; i<m; i++){
		for(int j=0; j<i; j++){
			R_ij = gsl_matrix_get(R, i, j);
			if(R_ij!=0) status=false;
		}
	}

	printf(status ? "True\n": "False\n");
	
	printf("checking Q^TQ=I by using gsl_BLAS:\n:");

	gsl_matrix* A_T = gsl_matrix_calloc(m,n);

	status_int = gsl_matrix_transpose_memcpy(A_T, A);

	gsl_matrix* C = gsl_matrix_calloc(m, m);

	status_int = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_T, A, 0.0, C);

	gsl_matrix* I = gsl_matrix_calloc(m, m);
	gsl_matrix_set_identity(I);

	printf("This is Q^TQ:\n");
	printm(C);

	gsl_matrix* leftside = gsl_matrix_calloc(n, m);

	status_int = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, leftside);

	int equal_final = gsl_matrix_equal(leftside, A_clone);

	printf("checking QR=A by using gsl_BLAS:\n:");
	
	printf("This is QR:\n");
	printm(leftside);
	
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(A_clone);

	gsl_matrix_free(A_T);
	gsl_matrix_free(C);
	gsl_matrix_free(I);
	gsl_matrix_free(leftside);
	
	return EXIT_SUCCESS;
}