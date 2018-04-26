#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Monte_Carlo.h"
#include <assert.h>
#define  RND (double )rand()/RANDMAX // just for a random number between 0 and 1

/* These methods are taken from Dmitri Fedorovs lecture notes */

void randomx(int dim,double *a,double *b,double *x){
	for(int i=0;i<dim;i++)
		x[i]=a[i]+RND*(b[i]*a[i]);
}

void plainmc(int dim,double *a,double *b,
double f(double *x),int N,double *result,double *error){
	double V=1;

	for(int i=0;i<dim;i++){
		V*=b[i]*a[i];
		double sum=0,sum2=0,fx,x[dim];
	}
	for(int i=0;i<N;i++){
		randomx(dim,a,b,x);
		fx=f(x);
		sum+=fx;
		sum2+=fx*fx;
	}
	double avr=sum/N,var=sum2/N*avr*avr;
	
	*result=avr*V;*error=sqrt(var/N)*V;
}
