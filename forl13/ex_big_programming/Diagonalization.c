#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <errno.h>
#define RND (double)rand()/RAND_MAX

double diagonalization(int n) {

	gsl_matrix *M = gsl_matrix_alloc(n, n);

	for(int i=0; i<n; i++) {
		for (int j=0; j<=i; j++) {

			double rand_number = RND;
			
			gsl_matrix_set(M, i, j, rand_number);
			gsl_matrix_set(M, j, i, rand_number);
		} 
	}



	/* for printing the matrix i case of emergency
	for(int i=0; i<n; i++) {
		printf("\n");
		for (int j=0; j<n; j++) {
			double a= gsl_matrix_get(M, i, j);
			printf("%g \t", a);
		} 
	}
	*/

	/* defining workspace and finding eigenvalues*/
	gsl_eigen_symm_workspace *work = gsl_eigen_symm_alloc(n); 

	int status; 
	gsl_vector *eval = gsl_vector_alloc(n);

	status = gsl_eigen_symm(M, eval, work);

	if(status != GSL_SUCCESS) {
		printf("Status = %s \n",  gsl_strerror (status));
	}

	/* print out the eigenvector 
	for(int i=0; i<n; i++) {
		printf("%g\n", gsl_vector_get(eval, i));
	}
	*/
	
	/* sorting and finding the largest */
	gsl_sort_vector(eval);
	double largest_eigenvalue = gsl_vector_get(eval, n-1);

	/* the largest 
	printf("-------------%g", largest_eigenvalue);
	*/

	/* frigør det hele */
	gsl_eigen_symm_free(work);
	gsl_matrix_free(M);
	gsl_vector_free(eval);

	return largest_eigenvalue;
}

int main() {

	const int max_n = 100;

	/* finding the biggest eigenvalue for several n
	   and directly printing it out. */
	for(int n=1; n<=max_n; n++) {

		double eigval = diagonalization(n);
		printf("%i \t %g \n", n, eigval);

	}
	return EXIT_SUCCESS;
}