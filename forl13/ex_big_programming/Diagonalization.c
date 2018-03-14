#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <errno.h>
#include <gsl/gsl_fit.h>
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

	double n_list[max_n];
	double eigval_list[max_n];

	/* finding the biggest eigenvalue for several n
	   and directly printing it out. */

	FILE *output  =fopen("output.txt", "w+");
	for(int n=1; n<=max_n; n++) {

		double eigval = diagonalization(n);
		fprintf(output, "%i \t %g \n", n, eigval);
		
		n_list[n-1]=(double) n;
		eigval_list[n-1]=eigval;

	}
	fclose(output);

	double c1, cov11, sumsq;
	double xstride=1, ystride=1;
	int status = gsl_fit_mul (n_list, xstride, \
		eigval_list, ystride, max_n, &c1, &cov11, &sumsq);
	
	if(status != GSL_SUCCESS) {
		printf("Status = %s \n",  gsl_strerror (status));
	}

	printf("this is c1=%g\n", c1);
	FILE *myfile= fopen("fit.txt","w+");
	for(int i=1; i<=100; i +=1) {
		fprintf(myfile, "%i\t %g\n", i, c1*i);
	}
	fclose(myfile);


	return EXIT_SUCCESS;
}