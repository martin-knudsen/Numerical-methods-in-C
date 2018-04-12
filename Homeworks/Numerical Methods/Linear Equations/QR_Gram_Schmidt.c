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

void qr_gs_inverse(const gsl_matrix* Q, const gsl_matrix* R, gsl_matrix* B){
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

void qr_gs_decomp_giddens(gsl_matrix* A, gsl_matrix* R) {
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

/* taken from Dmitri Fedorovs implementation.

*/
void givens_qr(gsl_matrix* A){ /*A<-Q,R*/
	// Iterating through the bottom part of matrix A first columns then rows
	for(int q=0;q<A->size2;q++)for(int p=q+1;p<A->size1;p++){
		// Finding out the angle that will zero out xq
		double theta=atan2(gsl_matrix_get(A,p,q),gsl_matrix_get(A,q,q));
		// iterating over above diagonal part of A
		for(int k=q;k<A->size2;k++){
			// getting the current value of xq, xp
			double xq=gsl_matrix_get(A,q,k),xp=gsl_matrix_get(A,p,k);
			// changing xq, xp according to the transformation and saving theta at A_pq
			gsl_matrix_set(A,q,k,xq*cos(theta)+xp*sin(theta));
			gsl_matrix_set(A,p,k,-xq*sin(theta)+xp*cos(theta));}
			gsl_matrix_set(A,p,q,theta);}}

// successively applying the givens transformation by using the found theta angles
void givens_qr_QTvec(gsl_matrix* QR,gsl_vector *v){/*v<-QˆTv*/
	// iterating columns of recently found QR 
	for(int q=0;q<QR->size2;q++)for(int p=q+1;p<QR->size1;p++){
		// saving the found theta value
		double theta=gsl_matrix_get(QR,p,q);
		//changing the appropriate entrance of v according to the giddens transformation
		double vq=gsl_vector_get(v,q),vp=gsl_vector_get(v,p);
		gsl_vector_set(v,q,vq*cos(theta)+vp*sin(theta));
		gsl_vector_set(v,p,-vq*sin(theta)+vp*cos(theta));}}

void givens_qr_solve(gsl_matrix* QR,gsl_vector* b){
	givens_qr_QTvec(QR,b);
	
	/* backsubstitution from Dmitri Fedorovs following (2.4)
		first time skipping the first for loop for starting point.
		then iterating down through yi each time first getting ci then
		subtracting the sum of Uik yik and dividing by Uii. 

	*/
	for(int i=b->size-1; i>=0; i--){
		double s = gsl_vector_get(b,i);
		for(int k=i+1; k<b->size; k++){
			s -= gsl_matrix_get(QR,i,k)*gsl_vector_get(b,k);
		}
		gsl_vector_set(b,i,s/gsl_matrix_get(QR,i,i));
	}
}

/* Also taken from Dmitri Fedorov. Basically using the same 
logic as my own Gram Schmidt implementation except now the 
system of equations is solved using givens
*/
void givens_qr_inverse(gsl_matrix* QR,gsl_matrix* B){
	gsl_matrix_set_identity(B);
	for(int i=0;i<QR->size2;i++){
		gsl_vector_view v=gsl_matrix_column(B,i);
	givens_qr_solve(QR,&v.vector);}}


int main() {
	
	const int n = 5;
	const int m = 4;

	gsl_matrix* A = gsl_matrix_alloc(n, m);
	gsl_matrix* A_clone = gsl_matrix_alloc(n, m);
	gsl_matrix* R = gsl_matrix_alloc(m, m);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			gsl_matrix_set(A, i, j, RND);
		}
	}

	gsl_matrix_memcpy(A_clone, A);
	qr_gs_decomp(A, R);

	printf("This is A:\n");
	printm(A_clone);
	printf("This is Q:\n");
	printm(A);
	printf("This is R:\n");
	printm(R);
	printf("Checking that R is upper triangular\n");
	printf("By checking wether left-down is zero \n");

	double R_ij; bool status=true;
	for(int i=1; i<m; i++){
		for(int j=0; j<i; j++){
			R_ij = gsl_matrix_get(R, i, j);
			if(R_ij!=0) status=false;
		}
	}

	printf(status ? "True\n": "False\n");
	
	printf("checking Q^TQ=I by using gsl_BLAS:\n:");

	gsl_matrix* A_T = gsl_matrix_alloc(m,n);

	gsl_matrix_transpose_memcpy(A_T, A);

	gsl_matrix* C = gsl_matrix_alloc(m, m);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_T, A, 0.0, C);

	gsl_matrix* I = gsl_matrix_alloc(m, m);
	gsl_matrix_set_identity(I);

	printf("This is Q^TQ:\n");
	printm(C);

	gsl_matrix* leftside = gsl_matrix_alloc(n, m);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, R, 0.0, leftside);


	printf("checking QR=A by using gsl_BLAS:\n:");
	
	printf("This is QR:\n");

	printm(leftside);

	printf("2. new matrix B:\n");
	int n_new  = 4;
	gsl_matrix* B = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix* B_clone = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix* R2 = gsl_matrix_alloc(n_new, n_new);

	for(int i=0; i<n_new; i++){
		for(int j=0; j<n_new; j++){
			double a = RND;
			gsl_matrix_set(B, i, j, a);
			gsl_matrix_set(B_clone, i, j, a);
		}
	}

	qr_gs_decomp(B, R2);

	printm(B_clone);
	printf("This is Q:\n");
	printm(B);
	printf("This is R:\n");
	printm(R2);

	gsl_vector* x = gsl_vector_alloc(n_new);
	gsl_vector* b = gsl_vector_alloc(n_new);
	gsl_vector* b_leftside = gsl_vector_alloc(n_new);
	
	printf("The vector b is:\n");
	for(int j=0; j<n_new; j++){
			gsl_vector_set(b, j, RND);
		}

	printv(b);

	qr_gs_solve(B, R2, b, x);
	
	printf("This is the solution x:\n");
	printv(x);

	gsl_blas_dgemv(CblasNoTrans, 1.0, B_clone, x, 0.0, b_leftside);

	printf("The vector Ax is:\n");
	printv(b_leftside);

	// part two
	printf("We'll use the same matrix for the inverse\n");
	printf("so A^-1:\n");
	gsl_matrix* D = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix* D_inverse = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix* D_test = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix_memcpy(D, B_clone);
	// remember B=Q, R2 is R
	qr_gs_inverse(B, R2, D_inverse);
	printm(D_inverse);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, D_inverse, 0.0, D_test);
	printf("AA^-1:\n");
	printm(D_test);

	//part three

	gsl_matrix* A_g = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix_memcpy(A_g, B_clone);
	gsl_matrix*  A_g_clone= gsl_matrix_alloc(n_new, n_new);
	gsl_matrix_memcpy(A_g_clone, B_clone);
	

	givens_qr(A_g);

	printf("3. givens stuff\n");
	printf("This is A_g:\n");
	printm(A_g_clone);
	printf("This is QR:\n");
	printm(A_g);
	printf("This is a new b:\n");

	gsl_vector* b_g = gsl_vector_alloc(n_new);
	gsl_vector* b_g_clone = gsl_vector_alloc(n_new);
	gsl_vector* A_gx = gsl_vector_alloc(n_new);
	for(int i=0; i<n_new; i++) {
		double a = RND;
		gsl_vector_set(b_g, i,a);
		gsl_vector_set(b_g_clone, i,a);
	}
	printv(b_g);


	givens_qr_solve(A_g, b_g);

	printf("this is the solution Ax=b using givens:\n");
	printv(b_g);

	printf("checking that this is correct by Ax:\n");
	gsl_blas_dgemv(CblasNoTrans, 1.0, A_g_clone, b_g, 0.0, A_gx);

	printv(A_gx);

	
	// part two
	printf("We'll use the same matrix for the inverse\n");
	printf("so A_G^-1:\n");
	gsl_matrix* D_g = gsl_matrix_alloc(n_new, n_new);
	gsl_matrix* D_g_test = gsl_matrix_alloc(n_new, n_new);
	
	// remember B=Q, R2 is R
	givens_qr_inverse(A_g,D_g);
	printm(D_g);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A_g_clone, D_g, 0.0, D_g_test);
	printf("AA^-1:\n");
	printm(D_g_test);
	
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(A_clone);
	gsl_matrix_free(B);
	gsl_matrix_free(R2);
	gsl_matrix_free(B_clone);
	gsl_matrix_free(D);
	gsl_matrix_free(D_inverse);
	gsl_matrix_free(D_test);
	gsl_matrix_free(A_g);
	gsl_matrix_free(A_g_clone);
	gsl_matrix_free(D_g);
	gsl_matrix_free(D_g_test);


	gsl_matrix_free(A_T);
	gsl_matrix_free(C);
	gsl_matrix_free(I);
	gsl_matrix_free(leftside);

	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_vector_free(b_leftside);
	gsl_vector_free(b_g);
	gsl_vector_free(b_g_clone);
	gsl_vector_free(A_gx);
	
	return EXIT_SUCCESS;
}