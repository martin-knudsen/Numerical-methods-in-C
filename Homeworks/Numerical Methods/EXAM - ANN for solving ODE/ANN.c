#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
typedef struct {int n; double (*g)(double); double (*gm)(double); gsl_vector* data;} ann;

ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double),double(*activation_function_der)(double)){
	ann* ann1=malloc(sizeof(ann));
	ann1->n=number_of_hidden_neurons;
	ann1->g=activation_function;
	ann2->gm=activation_function_der;
	ann1->data=gsl_vector_alloc(3*number_of_hidden_neurons);
	return ann1;
}

void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}

void ann_feed_forward(ann* network, double* F, double* dF){
		
	double argument,result_F, result_dF,sum_F=0,sum_dF=0, ai, bi, wi;
	for(int i=0; i<network->n;i++){
		
		ai = gsl_vector_get(network->data,i*3);
		bi = gsl_vector_get(network->data,i*3+1);
		wi = gsl_vector_get(network->data,i*3+2);
		argument = (x+ai)/bi;
		result_F = network->g(argument)*wi;
		result_dF = network->gm(argument)*wi/bi;
		sum_F += result_F;
		sum_dF += result_dF;
	}
	*F = sum_F;
	*dF = sum_dF;
}

void ann_train(ann* network, gsl_vector* xlist,double (*function)(double) ,double x0, double y0){
	int N=xlist->size;
	int n=network->n;
	double epsilon=1e-6, dx=1e-6, sum=0,Fpk,Fmpk,fk,xk,Fp0, Fmp0;

	double deviation_function (gsl_vector* p){
		gsl_vector_memcpy(network->data,p);
		sum=0;
		for(int k=0;k<N;k++){
			xk = gsl_vector_get(xlist,k);
			ann_feed_forward(network, xk,&Fpk,&Fmpk);
			fk = function(xk,Fpk)
			ann_feed_forward(network,x0,&Fp0,&Fmp0);
			sum += (Fmpk-f)*(Fmpk-f)+N*(Fp0-y0)*(Fp0-y0);
		}
		return sum;//N
	}
	gsl_vector* p=gsl_vector_alloc(3*n);
	gsl_vector_memcpy(p,network->data);
	newton_minimization_Broyden_num_gradient(deviation_function,p,dx,epsilon);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}
