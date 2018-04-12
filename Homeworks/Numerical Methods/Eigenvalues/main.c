#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f" //format of print "7 width, 3 digits after comma" 
#include "jacobi_cyclic.h"
#include "jacobi_eig_by_eig.h"
#include "jacobi_classic.h"

// inspired by Dmitri Fedorovs print implementation
void printm(gsl_matrix *A){
	// iterating over all rows and columns 
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) {
			printf(FMT,gsl_matrix_get(A,i,j));
		}
		printf("\n");}
}

void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}


int main() {
		
	const int n = 3;

	gsl_matrix* A = gsl_matrix_alloc(n, n);
	gsl_matrix* A_clone = gsl_matrix_alloc(n, n);
	gsl_matrix* V = gsl_matrix_alloc(n, n);
	gsl_matrix* D = gsl_matrix_alloc(n, n);
	gsl_matrix*	VTAV = gsl_matrix_alloc(n, n);
	gsl_vector* e = gsl_vector_alloc(n);

	gsl_matrix_set_identity(D);

	double a;
	for(int i=0; i<n; i++){
		for(int j=i; j<n; j++){
			a = RND;
			gsl_matrix_set(A, i, j, a);
			gsl_matrix_set(A, j, i, a);
		}
	}

	gsl_matrix_memcpy(A_clone, A);

	// start 
	int sweaps = jacobi_cyclic(A,e,V);

	for(int i=0; i<n; i++){
			gsl_matrix_set(D, i, i, gsl_vector_get(e, i));
	}

	printf("This is A:\n");
	printm(A_clone);
	printf("Using cyclic jacobi it has eigenvalue matrix D:\n");
	printm(D);
	printf("and eigenvectors:\n");
	printm(V);
	printf("It took only %i sweaps\n", sweaps);
	printf("Checking that V^TAV:\n");

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_clone, V, 0.0, A);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, V, A, 0.0, VTAV);
	printm(VTAV);

	// Part B 



	gsl_matrix_free(A);
	gsl_matrix_free(A_clone);
	gsl_matrix_free(V);
	gsl_matrix_free(VTAV);
	gsl_matrix_free(D);

	gsl_vector_free(e);
		
	return EXIT_SUCCESS;
}