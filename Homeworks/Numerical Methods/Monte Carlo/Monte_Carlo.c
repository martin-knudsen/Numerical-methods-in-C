#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Monte_Carlo.h"
#include <assert.h>
#define  RND (double )rand()/RANDMAX // just for a random number between 0 and 1

/* These methods are taken from Dmitri Fedorovs lecture notes 

	This function just */
void randomx(int dim,double *a,double *b,double *x){
	for(int i=0;i<dim;i++){
		x[i]=a[i]+RND*(b[i]-a[i]);
	}
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

/*
void plainmc2(double f(double *x),double  *a,double* b, N){
   var randomx = function(a,b) [a[i]+Math.random()*(b[i]-a[i]) for (i in a)];
   var volume=1; for(var i in a) volume*=b[i]-a[i];
   var sum=0,sum2=0;
   for(var i=0;i<N;i++){var fx=f(randomx(a,b)); sum += fx; sum2 += fx*fx}
   var mean = sum/N;                             // <f_i>
   var sigma = Math.sqrt(sum2/N - mean*mean);    // sigma² = <(f_i)²> - <f_i>²
   var SIGMA = sigma/Math.sqrt(N);               // SIGMA² = <Q²> - <Q>²
   return [mean*volume, SIGMA*volume];
}
*/