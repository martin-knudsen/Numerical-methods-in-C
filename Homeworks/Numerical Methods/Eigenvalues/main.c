#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <time.h>
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
#include "QR_ls.h"

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
		
	int n = 5;

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
	int number_rot_cyclic = 0;

	// start 
	int sweaps = jacobi_cyclic(A,e,V,&number_rot_cyclic);

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

	clock_t start, end;
	double cpu_time_used;

	int n_start=2, n_max=20, n_delta=1, m=1;
	n=(n_max-n_start)/n_delta+1;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* dy = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* COV = gsl_matrix_alloc(m, m);
	FILE* timedata=fopen("cyclic_time.txt","w+");
	int i=0;
	for(n=n_start;n<=n_max;n+=n_delta) {
		gsl_matrix* A_time = gsl_matrix_alloc(n, n);
		gsl_matrix* V_time = gsl_matrix_alloc(n, n);
		gsl_vector* e_time = gsl_vector_alloc(n);
		for(int i=0; i<n; i++){
			for(int j=i; j<n; j++){
				a = RND;
				gsl_matrix_set(A_time, i, j, a);
				gsl_matrix_set(A_time, j, i, a);
			}
		}
		double number_rot_cyclic2=0;
		start = clock();
		sweaps = jacobi_cyclic(A_time,e_time,V_time,&number_rot_cyclic2);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		fprintf(timedata,"%i\t%g\n",n,cpu_time_used);
		gsl_vector_set(x,i,n);
		gsl_vector_set(y,i,cpu_time_used);
		gsl_vector_set(dy,i,1.0);
		gsl_matrix_free(A_time);
		gsl_matrix_free(V_time);
		gsl_vector_free(e_time);
		i++;
	}
	fclose(timedata);

	double funs2(int i, double x){
   		switch(i){
   		case 0: return x*x*x; break;
	   	default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   		}
	}
	
	least_squares(x,y,dy,m,funs2,c,COV);
	FILE* errorfit = fopen("cyclic_time_fit.txt","w+");
	double delta_N=0.2, max_N=n_max; 
	for(double N=1;N<max_N;N+=delta_N){
		fprintf(errorfit,"%g\t%g\n",N,gsl_vector_get(c,0)*N*N*N);
	}
	fclose(errorfit);
	printf("\nCheck that error behaves as O(n^3) by Least squares fit\n");
	printf("of the cyclic jacobi algorithm from n=%i to n=%i\n",n_start,n_max);	
	printf("Fit coefficient found:%g\n",gsl_vector_get(c,0));
	printf("As one can see on the plot.svg the time follows very nicely\n"
		"the O(n^3) prediction. There is of course uncertainty from the\n"
		"prediction because of the randomness component in initializing matrix A\n\n");
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_matrix_free(COV);
	gsl_vector_free(c);
	

	// Part B 
	int n_eigval = 1; int number_rot_eig_by_eig = 0;
	printf("1. Using the same matrix A as before for comparison\n");
	printf("The resulting %i lowest eigenvalues are\n", n_eigval);
	gsl_matrix_set_identity(V);
	gsl_vector_set_zero(e);

	gsl_matrix_memcpy(A, A_clone);
	sweaps = jacobi_eig_by_eig(A,e,V, n_eigval, &number_rot_eig_by_eig);

	for(int i=0; i<n_eigval; i++) {
		printf(FMT, gsl_vector_get(e, i));
		printf("\n");
	}
	printf("The reason for this being the lowest is that this \n");
	printf("Eigenvalue is calculated by the most subtractions of all\n");
	printf("In the algorithm according to (3.9.4) \n");

	n_eigval = 2; int number_rot_waste = 0;
	printf("The resulting %i lowest eigenvalues are\n", n_eigval);
	gsl_matrix_set_identity(V);
	gsl_vector_set_zero(e);

	gsl_matrix_memcpy(A, A_clone);
	sweaps = jacobi_eig_by_eig(A,e,V, n_eigval, &number_rot_waste);

	for(int i=0; i<n_eigval; i++) {
		printf(FMT, gsl_vector_get(e, i));
		printf("\n");
	}
	printf("The algorithm works by constantly making the off-diagonal elements \n");
	printf("smaller. when we have already zeroed a row\n");
	printf("there is no reason to make it smaller\n");
	printf("Just as before we have the second most subtraction of any diagonal\n");
	printf("element\n");

	n_eigval = 3;
	printf("The resulting %i lowest eigenvalues are\n", n_eigval);
	gsl_matrix_set_identity(V);
	gsl_vector_set_zero(e);
	number_rot_waste = 0;
	gsl_matrix_memcpy(A, A_clone);
	sweaps = jacobi_eig_by_eig(A,e,V, n_eigval, &number_rot_waste);

	for(int i=0; i<n_eigval; i++) {
		printf(FMT, gsl_vector_get(e, i));
		printf("\n");
	}

	printf("2. The way to find the biggest eigenvalues is just \n");
	printf("to add pi/2 to phi\n");

	printf("normal run of the algorithm takes %i rotations\n", number_rot_cyclic);
	printf("Finding only the first eigenvalue using eig by eig takes %i rotations\n", number_rot_eig_by_eig);
	printf("Finding all eigenvalues using eig by eig takes: %i\n", number_rot_waste);
	printf("So it doesn't take less rotations, but maybe less time because of a saved for-loop\n");

	printf("C. The classic algorithm gives this result:\n");
	number_rot_waste = 0;
	gsl_matrix_memcpy(A, A_clone);
	gsl_matrix_set_identity(V);
	gsl_vector_set_zero(e);
	sweaps = jacobi_classic(A,e,V, &number_rot_waste);
	printv(e);
	printf("This many rotations %i\n",number_rot_waste);
	printf("So it is faster than the cyclical as expected. However this comes at the cost of\n"
		"having to update the list of indices. This however is worth it since it takes O(n) time\n"
		"and is therefore not the bottleneck of the algorithm (atleast if n is big)\n");
	gsl_matrix_free(A);
	gsl_matrix_free(A_clone);
	gsl_matrix_free(V);
	gsl_matrix_free(VTAV);
	gsl_matrix_free(D);

	gsl_vector_free(e);
		
	return EXIT_SUCCESS;
}