#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Monte_Carlo.h"
#include <assert.h>
#define  RND (double )rand()/RAND_MAX // just for a random number between 0 and 1

/* These methods are taken from Dmitri Fedorovs lecture notes 

	This function just calcuates a random number within the integration limits
	for each integration variable, meaning each dimension of array x */
void randomx(int dim,double *a,double *b,double *x){
	// for each dimension rescale the RND value to be in the interval of the limits for that variable
	// and assigns them to the x array
	for(int i=0;i<dim;i++){
		x[i]=a[i]+RND*(b[i]-a[i]);
	}
}

/* This is the actual plain Monte Carlo integrator. 
It takes in the dimension of the problem, the limits of all the integration variables
as arrays a and b, the integration function, taking an array x of variables.
The number of random points N, the result and the error*/
void plainmc(int dim,double *a,double *b,
double f(double *x),int N,double *result,double *error){
	double V=1;
	// Here we set the multidiminsional volume as a multidimensional "rectangle"
	// containing the integration volume Omega. Just by multiplying with sides of the 
	// rectangle as big as the integration limits for each integration variable.
	for(int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}
	double sum=0,sum2=0,fx,x[dim];
	
	// loop for finding random multi-Dimensional point, getting the 
	// function value at that point, and calculating the sum of the 
	// function values and the function values squared for finding the
	// expectation value of both.
	for(int i=0;i<N;i++){
		randomx(dim,a,b,x);
		fx=f(x);
		sum+=fx;
		sum2+=fx*fx;
	}
	// Average <f> by dividing with number of points
	double avr=sum/N;
	//finding the variance, sigma^2, by using the standard formular (9.4) 
	double var=sum2/N-avr*avr;
	// estimate of the integral is just (9.2) 
	*result=avr*V;
	// the error is just the squareroot of the variance
	*error=sqrt(var/N)*V;
}

