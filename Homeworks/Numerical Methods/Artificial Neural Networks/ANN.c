#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "Newton.h"
typedef struct {int n; double (*f)(double); gsl_vector* data;} ann;
typedef struct {int n; double (*f)(double); gsl_vector* data;} ann2D;

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
	gsl_vector_memcpy(p,network->data);
	newton_minimization_Broyden_num_gradient(func,p,dx,epsilon);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}

ann2D* ann2D_alloc(int number_of_hidden_neurons, double(*activation_function)(double)){
	ann2D* ann1=malloc(sizeof(ann2D));
	ann1->n=number_of_hidden_neurons;
	ann1->f=activation_function;
	ann1->data=gsl_vector_alloc(5*number_of_hidden_neurons);
	return ann1;
}

void ann2D_free(ann2D* network){
	gsl_vector_free(network->data);
	free(network);
}

double ann2D_feed_forward(ann2D* network, double x1, double x2){
	
	double argument1,argument2,result,sum=0,ai1,bi1,wi1,ai2,bi2,wi2;
	for(int i=0; i<network->n;i++){
		
		ai1=gsl_vector_get(network->data,i*5);
		//printf("%g\n",ai1);
		bi1=gsl_vector_get(network->data,i*5+1);
		//printf("%g\n",bi1);
		argument1=(x1+ai1)/bi1;
		//printf("%g\n",argument1);
		ai2=gsl_vector_get(network->data,i*5+2);
		//printf("%g\n",ai2);
		bi2=gsl_vector_get(network->data,i*5+3);
		//printf("%g\n",bi2);
		wi2=gsl_vector_get(network->data,i*5+4);
		//printf("%g\n",wi2);
		argument2=(x2+ai2)/bi2;
		//printf("%g\n",argument2);
		result=(network->f(argument1)*network->f(argument2))*wi2;
		//printf("%g\n",result);
		sum+=result;
		
	}
	return sum;

}
void ann2D_train(ann2D* network, gsl_matrix* xlist, gsl_vector* ylist){
	int N=ylist->size;
	int n=network->n;
	double epsilon=1e-6, dx=1e-3, sum=0,Fpk,xk1,xk2,yk;

	double func (gsl_vector* p){
		gsl_vector_memcpy(network->data,p);
		sum=0;
		for(int k=0;k<N;k++){
			xk1=gsl_matrix_get(xlist,k,0);
			xk2=gsl_matrix_get(xlist,k,1);
			yk=gsl_vector_get(ylist,k);
			Fpk=ann2D_feed_forward(network, xk1,xk2);
			sum+=(Fpk-yk)*(Fpk-yk);
		}
		//printf("%g\n",sum);
	
		return sum/N;
	}
	gsl_vector* p=gsl_vector_alloc(5*n);
	gsl_vector_memcpy(p,network->data);
	//printv(p);

	newton_minimization_Broyden_num_gradient(func,p,dx,epsilon);
	gsl_vector_memcpy(network->data,p);
	//printv(p);
	gsl_vector_free(p);
}

/*

double ann2D_feed_forward(ann2D* network, double x1, double x2){
	
	double argument1,argument2,result,sum=0,ai1,bi1,wi1,ai2,bi2,wi2;
	for(int i=0; i<network->n;i++){
		
		ai1=gsl_vector_get(network->data,i*5);
		//printf("%g\n",ai1);
		bi1=gsl_vector_get(network->data,i*5+1);
		//printf("%g\n",bi1);
		argument1=(x1+ai1)/bi1;
		//printf("%g\n",argument1);
		ai2=gsl_vector_get(network->data,i*5+2);
		//printf("%g\n",ai2);
		bi2=gsl_vector_get(network->data,i*5+3);
		//printf("%g\n",bi2);
		wi2=gsl_vector_get(network->data,i*5+4);
		//printf("%g\n",wi2);
		argument2=(x2+ai2)/bi2;
		//printf("%g\n",argument2);
		result=network->f(argument1,argument2)*wi2;
		//printf("%g\n",result);
		sum+=result;
		
	}
	return sum;

}
void ann2D_train(ann2D* network, gsl_matrix* xlist, gsl_vector* ylist){
	int N=ylist->size;
	int n=network->n;
	double epsilon=1e-6, dx=1e-6, sum=0,Fpk,xk1,xk2,yk;

	double func (gsl_vector* p){
		gsl_vector_memcpy(network->data,p);
		sum=0;
		for(int k=0;k<N;k++){
			xk1=gsl_matrix_get(xlist,k,0);
			xk2=gsl_matrix_get(xlist,k,1);
			yk=gsl_vector_get(ylist,k);
			Fpk=ann2D_feed_forward(network, xk1,xk2);
			sum+=(Fpk-yk)*(Fpk-yk);
		}
		//printf("%g\n",sum);
		return sum/N;
	}
	gsl_vector* p=gsl_vector_alloc(5*n);
	gsl_vector_memcpy(p,network->data);
	//printv(p);

	qnewton(func,p,epsilon);
	gsl_vector_memcpy(network->data,p);
	printv(p);
	gsl_vector_free(p);
}

*/
/*
ann* ann_digits_alloc(int number_of_hidden_neurons, double(*activation_function)(double)){
	ann* ann1=malloc(sizeof(ann));
	ann1->n=number_of_hidden_neurons;
	ann1->f=activation_function;
	ann1->data=gsl_vector_alloc(3*number_of_hidden_neurons*35);
	return ann1;
}

void ann_feed_forward(ann* network, gsl_vector* x,gsl_vector* probs){
		

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
	gsl_vector_memcpy(p,network->data);
	newton_minimization_Broyden_num_gradient(func,p,dx,epsilon);
	gsl_vector_memcpy(network->data,p);
	gsl_vector_free(p);
}

*/