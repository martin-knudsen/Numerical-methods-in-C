#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "QR.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f"

void jacobian_num(void f(gsl_vector* x, gsl_vector* fx),
				  gsl_vector* x, gsl_matrix* J, double dx){
int n = J->size1; 

gsl_vector* fx = gsl_vector_alloc(n);

double Jik, fi_1, fi_2, xk;
for(int i=0; i<n; i++){
	for(int k=0; k<n; k++){
		f(x,fx);
		fi_1=gsl_vector_get(fx, i);
		xk = gsl_vector_get(x,k);
		gsl_vector_set(x,k,xk+dx);
		f(x,fx);
		gsl_vector_set(x,k,xk-dx);
		fi_2=gsl_vector_get(fx,i);
		Jik=(fi_2-fi_1)/dx;
		gsl_matrix_set(J,i,k,Jik);
	}
}

gsl_vector_free(fx);

}

void newton_num(
	void f(gsl_vector* x, gsl_vector* fx),
	gsl_vector* x,
	double dx,
	double epsilon
){

	int n = x->size;
	gsl_vector* x_plus = gsl_vector_alloc(n);
	gsl_vector* deltax = gsl_vector_alloc(n);
	gsl_vector* fx_minus = gsl_vector_alloc(n);
	gsl_vector* fx_plus = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_matrix* J = gsl_matrix_alloc(n,n);

	

	double fx_norm=INFINITY, deltax_norm;
	double lambda, fx_plus_norm; 
	while(fx_norm>epsilon){
		f(x,fx);
		gsl_vector_memcpy(fx_minus, fx);
		gsl_vector_scale(fx_minus, -1.0);
		gsl_vector_memcpy(x_plus,x);

		jacobian_num(f,x,J,dx);
		qr_solve_system(J, fx_minus, deltax);

		lambda=1.0;
		gsl_vector_scale(x_plus,lambda);
		gsl_vector_add(x_plus, x);
		f(x_plus,fx_plus);
		fx_norm=gsl_blas_dnrm2(fx);
		fx_plus_norm=gsl_blas_dnrm2(fx_plus);
		deltax_norm=gsl_blas_dnrm2(deltax);

		while(fx_plus_norm>(1-lambda/2)*fx_norm 
			  && deltax_norm*lambda>dx){
			gsl_vector_scale(deltax,lambda);
			gsl_vector_add(x, deltax);
		}
	}

	gsl_vector_free(x_plus);
	gsl_vector_free(deltax);
	gsl_vector_free(fx);
	gsl_vector_free(fx_minus);
	gsl_vector_free(fx_plus);
	gsl_matrix_free(J);
}