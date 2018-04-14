#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "jacobi_cyclic.h"

void singular_decomp(gsl_matrix* A, gsl_matrix* U, 
					 gsl_matrix* S, gsl_matrix* V) {
	int n = A->size1, m = A->size2;

	gsl_matrix* ATA = gsl_matrix_alloc(m,m);
	gsl_matrix* D = gsl_matrix_alloc(m,m);
	gsl_matrix* D_minus = gsl_matrix_alloc(m,m);
	gsl_matrix* AV = gsl_matrix_alloc(n,m);
	gsl_matrix* US = gsl_matrix_alloc(n,m);
	gsl_matrix* USVT = gsl_matrix_alloc(n,m);
	gsl_vector* e = gsl_vector_alloc(m);
	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, ATA);
	
	int number_rot=0;
	int sweeps =jacobi_cyclic(ATA,e,V,&number_rot); 
	
	double ei;
	for(int i=0; i<m; i++){
		ei = gsl_vector_get(e, i);
		gsl_matrix_set(D, i, i, ei);
		gsl_matrix_set(S, i, i, sqrt(ei));
		gsl_matrix_set(D_minus, i, i, 1.0/sqrt(ei));

	}
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, V, 0.0, AV);
	
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AV, D_minus, 0.0, U);
	
	// checking that this implementation is correct. 
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, S, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, V, 0.0, USVT);
	

	gsl_matrix_free(ATA);
	gsl_matrix_free(D);
	gsl_matrix_free(D_minus);
	gsl_matrix_free(AV);
	gsl_matrix_free(US);
	gsl_matrix_free(USVT);
	gsl_vector_free(e);
}

void least_squares_problem_singular(gsl_matrix* A, gsl_vector* x, 
									gsl_vector* b, gsl_matrix* COV) {
	int n = A->size1, m = A->size2;
	gsl_matrix* U = gsl_matrix_alloc(n,m);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	gsl_matrix* S_clone = gsl_matrix_alloc(m,m);
	gsl_matrix* S_minus2 = gsl_matrix_alloc(m,m);
	gsl_matrix* VS_minus2 = gsl_matrix_alloc(m,m);
	gsl_matrix* R = gsl_matrix_alloc(m,m);
	gsl_matrix* V = gsl_matrix_alloc(m,m);
	gsl_vector* UTb = gsl_vector_alloc(m);
	gsl_vector* y = gsl_vector_alloc(m);

	
	singular_decomp(A,U,S,V);
	
	gsl_matrix_memcpy(S_clone, S);

	// solving Sy=U^Tb
	gsl_blas_dgemv(CblasTrans,1.0,U,b,0.0,UTb);
	qr_gs_decomp(S_clone,R);
	// now S_clone is Q

	qr_gs_solve(S_clone,R,UTb,y);

	// have found least square solution x
	gsl_blas_dgemv(CblasNoTrans,1.0,V,y,0.0,x);

	// now will find COV
	// first set S^-2
	double Sii;
	for(int i=0;i<m;i++){
		Sii=gsl_matrix_get(S, i, i);
		gsl_matrix_set(S_minus2, i, i, 1.0/(Sii*Sii));
	}

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,S_minus2,0.0,VS_minus2);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,VS_minus2,V,0.0,COV);

	gsl_matrix_free(U);
	gsl_matrix_free(S);
	gsl_matrix_free(S_clone);
	gsl_matrix_free(R);
	gsl_matrix_free(V);
	gsl_matrix_free(S_minus2);
	gsl_matrix_free(VS_minus2);
	gsl_vector_free(UTb);
	gsl_vector_free(y);

}

void least_squares_singular(gsl_vector* x, gsl_vector* y, gsl_vector* dy, int m,
				 double funs(int i, double x), gsl_vector* c, gsl_matrix* COV){
	int n = x->size;
	gsl_matrix* A = gsl_matrix_alloc(n, m);
	gsl_matrix* A_clone = gsl_matrix_alloc(n, m);
	gsl_vector* b = gsl_vector_alloc(n);
	
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
		}
	}
	
	least_squares_problem_singular(A,c,b,COV);

	gsl_matrix_free(A);
	gsl_matrix_free(A_clone);
	gsl_vector_free(b);
}


	