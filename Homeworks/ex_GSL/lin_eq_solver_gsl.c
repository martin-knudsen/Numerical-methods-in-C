#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int main() {

	gsl_vector* y = gsl_vector_alloc(3);
	
	gsl_vector_set(y, 0, 6.23);
	gsl_vector_set(y, 1, 5.37);
	gsl_vector_set(y, 2, 2.29);

	gsl_matrix* A = gsl_matrix_alloc(3, 3);

	gsl_matrix_set(A, 0, 0, 6.13);
	gsl_matrix_set(A, 0, 1, -2.90);
	gsl_matrix_set(A, 0, 2, 5.86);
	gsl_matrix_set(A, 1, 0, 8.08);
	gsl_matrix_set(A, 1, 1, -6.31);
	gsl_matrix_set(A, 1, 2, -3.89);
	gsl_matrix_set(A, 2, 0, -4.36);
	gsl_matrix_set(A, 2, 1, 1.00);
	gsl_matrix_set(A, 2, 2, 0.19);

	gsl_vector* x = gsl_vector_alloc(3);

	gsl_linalg_HH_solve(A, y, x);

	double x1 = gsl_vector_get(x, 0);
	double x2 = gsl_vector_get(x, 1);
	double x3 = gsl_vector_get(x, 2); 

	printf("Trying HH method\n");
	printf("x = [%g %g %g]\n", x1, x2, x3);

	gsl_matrix* A_copy = gsl_matrix_alloc(3, 3);
	gsl_matrix_set(A_copy, 0, 0, 6.13);
	gsl_matrix_set(A_copy, 0, 1, -2.90);
	gsl_matrix_set(A_copy, 0, 2, 5.86);
	gsl_matrix_set(A_copy, 1, 0, 8.08);
	gsl_matrix_set(A_copy, 1, 1, -6.31);
	gsl_matrix_set(A_copy, 1, 2, -3.89);
	gsl_matrix_set(A_copy, 2, 0, -4.36);
	gsl_matrix_set(A_copy, 2, 1, 1.00);
	gsl_matrix_set(A_copy, 2, 2, 0.19);


	gsl_vector* y_test = gsl_vector_alloc(3);

	gsl_blas_dgemv(CblasNoTrans, 1.0, A_copy, x, 0, y_test);
	printf("Testing wether solution is correct.\n");
	
	double y1_test = gsl_vector_get(y_test, 0);
	double y2_test = gsl_vector_get(y_test, 1);
	double y3_test = gsl_vector_get(y_test, 2);

	double y1 = gsl_vector_get(y, 0);
	double y2 = gsl_vector_get(y, 1);
	double y3 = gsl_vector_get(y, 2);

	printf("[%g, %g, %g] = [%g %g %g]\n", \
		y1, y2, y3, y1_test, y2_test, y3_test);

	if(y1==y1_test && y2==y2_test && y3==y3_test) printf("Success!!\n");


	/* The optional part: make hilbert matrix*/

	gsl_matrix* H = gsl_matrix_alloc(4, 4);

	for(int i=0; i<4; i++) {
		for(int j=0; j<4; j++) {
			gsl_matrix_set(H, i, j, 1/(i+j+1));
		}
	}
	
	return 0;
}