#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "QR.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f"

void newton_minimization(
	double f(gsl_vector* x, gsl_vector* df, gsl_matrix* H),
	gsl_vector* x,
	double dx,
	double epsilon
){

	int n = x->size;
	gsl_vector* x_plus = gsl_vector_alloc(n);
	gsl_vector* deltax = gsl_vector_alloc(n);
	gsl_vector* df_minus = gsl_vector_alloc(n);
	gsl_vector* df_plus = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_matrix* H = gsl_matrix_alloc(n,n);
	
	double deltax_norm=INFINITY, df_norm= INFINITY, fx;
	double lambda, DxTdf, fx_plus, alfa=0.51; 
	int steps=0;
	while(df_norm>epsilon && deltax_norm>dx){
		steps++;
		// find f(x)
		fx = f(x,df,H);
		// convert to -f(x) for matrix stuff
		gsl_vector_memcpy(df_minus, df);
		gsl_vector_scale(df_minus, -1.0);
		
		qr_solve_system(H, df_minus, deltax);
		

		// lambda start value
		lambda=1.0;
		
		// make xplus (x+lambda*deltax);
		gsl_vector_memcpy(x_plus,deltax);
		gsl_vector_scale(x_plus,lambda);
		gsl_vector_add(x_plus, x);
		// now using that x_plus to create a new df_plus
		fx_plus=f(x_plus,df_plus,H);
		
		gsl_blas_ddot(deltax,df,&DxTdf);
		deltax_norm=gsl_blas_dnrm2(deltax);
		while(fx_plus>fx+alfa*lambda*DxTdf && lambda>1.0/64){
			lambda/=2;

			gsl_vector_memcpy(x_plus,deltax);
			gsl_vector_scale(x_plus,lambda);
			gsl_vector_add(x_plus, x);
			// now using that x_plus to create a new df_plus
			fx_plus=f(x_plus,df_plus,H);
		
			
			
		}
		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(x, deltax);
		fx=f(x, df,H);
		df_norm=gsl_blas_dnrm2(df);
	}
	printf("Number of steps taken is: %i\n", steps);

	gsl_vector_free(x_plus);
	gsl_vector_free(deltax);
	gsl_vector_free(df);
	gsl_vector_free(df_minus);
	gsl_vector_free(df_plus);
	gsl_matrix_free(H);
}