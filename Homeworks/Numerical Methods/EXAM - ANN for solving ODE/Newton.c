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

void gradient_num(double f(gsl_vector* x),
				  gsl_vector* x, gsl_vector* df, double dx){
int n = x->size; 

double dfi, fi_1, fi_2, xi;
for(int i=0; i<n; i++){
		fi_1=f(x);
		xi = gsl_vector_get(x,i);
		gsl_vector_set(x,i,xi+dx);
		fi_2=f(x);
		gsl_vector_set(x,i,xi-dx);
		dfi=(fi_2-fi_1)/dx;
		gsl_vector_set(df,i,dfi);
	}
}

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

void newton_minimization_Broyden(
	double f(gsl_vector* x, gsl_vector* df),
	gsl_vector* x,
	double dx,
	double epsilon
){

	int n = x->size;
	gsl_vector* x_plus = gsl_vector_alloc(n);
	gsl_vector* deltax = gsl_vector_alloc(n);
	gsl_vector* df_plus = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_matrix* H_inv = gsl_matrix_alloc(n,n);
	gsl_vector* H_inv_y = gsl_vector_alloc(n);
	gsl_vector* delta_H_inv_y = gsl_vector_alloc(n);
	gsl_matrix* d_H_inv = gsl_matrix_alloc(n,n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_matrix_set_identity(H_inv);

	double deltax_norm=INFINITY, df_norm= INFINITY, fx;
	double lambda, DxTdf, fx_plus, alfa=0.51, H_inv_update; 
	double update_num, update_denom;
	int steps=0;
	while(df_norm>epsilon && deltax_norm>dx){
		steps++;
		// find f(x)
		fx = f(x,df);
		gsl_vector_memcpy(x_plus,x);
		gsl_vector_add(x_plus,deltax);
		fx_plus = f(x_plus,df_plus);
		
		// To avoid updating the hessian in the first run. 
		if(steps>1){
		
		gsl_vector_memcpy(y,df_plus);
		gsl_vector_sub(y,df);
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,y,0,H_inv_y);
		gsl_vector_memcpy(delta_H_inv_y,deltax);
		gsl_vector_sub(delta_H_inv_y,H_inv_y);
		gsl_blas_ddot(delta_H_inv_y,deltax,&update_num);
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,deltax,0,H_inv_y);
		gsl_blas_ddot(y,H_inv_y,&update_denom);
		H_inv_update = update_num/update_denom;
		gsl_matrix_memcpy(d_H_inv, H_inv);
		gsl_matrix_scale(d_H_inv,H_inv_update);
		gsl_matrix_add(H_inv,d_H_inv);
		}
		// now we use all that to find deltax with the new H
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,df,0,deltax);
		gsl_vector_scale(deltax,-1.0);
		// lambda start value
		lambda=1.0;
		
		// make xplus (x+lambda*deltax);
		gsl_vector_memcpy(x_plus,deltax);
		gsl_vector_scale(x_plus,lambda);
		gsl_vector_add(x_plus, x);
		// now using that x_plus to create a new df_plus
		fx_plus=f(x_plus,df_plus);
		
		gsl_blas_ddot(deltax,df,&DxTdf);
		deltax_norm=gsl_blas_dnrm2(deltax);
		while(fx_plus>fx+alfa*lambda*DxTdf && lambda>1.0/64){
			lambda/=2;
			gsl_vector_memcpy(x_plus,deltax);
			gsl_vector_scale(x_plus,lambda);
			gsl_vector_add(x_plus, x);
			// now using that x_plus to create a new df_plus
			fx_plus=f(x_plus,df_plus);
			
		}
		
		// reset to unity if the update is negative.
		if(gsl_blas_dnrm2(df)<gsl_blas_dnrm2(df_plus)) {
			gsl_matrix_set_identity(H_inv);

		}
		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(x, deltax);
		fx=f(x, df);
		df_norm=gsl_blas_dnrm2(df);
	}
	printf("Number of steps taken is: %i\n", steps);
	gsl_vector_free(x_plus);
	gsl_vector_free(deltax);
	gsl_vector_free(df);
	gsl_vector_free(df_plus);
	gsl_matrix_free(H_inv);
	gsl_vector_free(y);
	gsl_vector_free(delta_H_inv_y);
	gsl_matrix_free(d_H_inv);
	gsl_vector_free(H_inv_y);

}

void newton_minimization_Broyden_num_gradient(
	double f(gsl_vector* x),
	gsl_vector* x,
	double dx,
	double epsilon
){
	int n = x->size;
	gsl_vector* x_plus = gsl_vector_alloc(n);
	gsl_vector* deltax = gsl_vector_alloc(n);
	gsl_vector* df_plus = gsl_vector_alloc(n);
	gsl_vector* df = gsl_vector_alloc(n);
	gsl_matrix* H_inv = gsl_matrix_alloc(n,n);
	gsl_vector* H_inv_y = gsl_vector_alloc(n);
	gsl_vector* delta_H_inv_y = gsl_vector_alloc(n);
	gsl_matrix* d_H_inv = gsl_matrix_alloc(n,n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_matrix_set_identity(H_inv);

	double deltax_norm=INFINITY, df_norm= INFINITY, fx;
	double lambda, DxTdf, fx_plus, alfa=0.51, H_inv_update; 
	double update_num, update_denom;
	int steps=0;
	while(df_norm>epsilon && deltax_norm>dx){
		steps++;
		// find f(x)
		fx = f(x);
		gradient_num(f,x,df,dx);
		gsl_vector_memcpy(x_plus,x);
		gsl_vector_add(x_plus,deltax);
		fx_plus = f(x_plus);
		gradient_num(f,x_plus,df_plus,dx);
		// To avoid updating the hessian in the first run. 
		if(steps>1){
		
		gsl_vector_memcpy(y,df_plus);
		gsl_vector_sub(y,df);
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,y,0,H_inv_y);
		gsl_vector_memcpy(delta_H_inv_y,deltax);
		gsl_vector_sub(delta_H_inv_y,H_inv_y);
		gsl_blas_ddot(delta_H_inv_y,deltax,&update_num);
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,deltax,0,H_inv_y);
		gsl_blas_ddot(y,H_inv_y,&update_denom);
		H_inv_update = update_num/update_denom;
		gsl_matrix_memcpy(d_H_inv, H_inv);
		gsl_matrix_scale(d_H_inv,H_inv_update);
		gsl_matrix_add(H_inv,d_H_inv);
		}
		// now we use all that to find deltax with the new H
		gsl_blas_dgemv(CblasNoTrans,1,H_inv,df,0,deltax);
		gsl_vector_scale(deltax,-1.0);
		// lambda start value
		lambda=1.0;
		
		// make xplus (x+lambda*deltax);
		gsl_vector_memcpy(x_plus,deltax);
		gsl_vector_scale(x_plus,lambda);
		gsl_vector_add(x_plus, x);
		// now using that x_plus to create a new df_plus
		fx_plus=f(x_plus);
		gradient_num(f,x_plus,df_plus,dx);
		
		gsl_blas_ddot(deltax,df,&DxTdf);
		deltax_norm=gsl_blas_dnrm2(deltax);
		while(fx_plus>fx+alfa*lambda*DxTdf && lambda>1.0/64){
			lambda/=2;
			gsl_vector_memcpy(x_plus,deltax);
			gsl_vector_scale(x_plus,lambda);
			gsl_vector_add(x_plus, x);
			// now using that x_plus to create a new df_plus
			fx_plus=f(x_plus);
			gradient_num(f,x_plus,df_plus,dx);
		}
		
		// reset to unity if the update is negative.
		if(gsl_blas_dnrm2(df)<gsl_blas_dnrm2(df_plus)) {
			gsl_matrix_set_identity(H_inv);

		}
		gsl_vector_scale(deltax,lambda);
		gsl_vector_add(x, deltax);
		fx=f(x);
		gradient_num(f,x,df,dx);
		df_norm=gsl_blas_dnrm2(df);
	}
	
	gsl_vector_free(x_plus);
	gsl_vector_free(deltax);
	gsl_vector_free(df);
	gsl_vector_free(df_plus);
	gsl_matrix_free(H_inv);
	gsl_vector_free(y);
	gsl_vector_free(delta_H_inv_y);
	gsl_matrix_free(d_H_inv);
	gsl_vector_free(H_inv_y);

}