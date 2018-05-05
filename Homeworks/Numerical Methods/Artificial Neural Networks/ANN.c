#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;


ann* ann_alloc(int number_of_hidden_neurons, double(*activation_function)(double)){
	ann* ann1=malloc(sizeof(ann));
	ann1->n=number_of_hidden_neurons;
	ann1->f=activation_function;
	ann1->data=gsl_vector_alloc(3*number_of_hidden_neurons);
	return ann1;
}

void ann_free(ann* network){
	gsl_vector_free(network->data);
	free(network);
}

double ann_feed_forward(ann* network, double x){
		

	double argument,result,sum=0, ai, bi, wi;
	for(int i=0; i<network->n;i++){
		
		ai=gsl_vector_get(network->data,i*3);
		bi=gsl_vector_get(network->data,i*3+1);
		wi=gsl_vector_get(network->data,i*3+2);
		argument=(x+ai)/bi;
		result=network->f(argument)*wi;
		sum+=result;
		
	}
	return sum;

}
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist){
	int N=xlist->size;
	int n=network->n;
	double epsilon=1e-6, dx=1e-6, sum=0,Fpk,xk,yk;

	double func (gsl_vector* p){
		gsl_vector_memcpy(network->data,p);
		sum=0;
		for(int k=0;k<N;k++){
			xk=gsl_vector_get(xlist,k);
			yk=gsl_vector_get(ylist,k);
			Fpk=ann_feed_forward(network, xk);
			sum+=(Fpk-yk)*(Fpk-yk);
		}
		return sum/N;
	}
	gsl_vector* p=gsl_vector_alloc(3*n);
	printf("1\n");
	gsl_vector_memcpy(p,network->data);printf("1\n");
	newton_minimization_Broyden_num_gradient(func,p,dx,epsilon);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}
