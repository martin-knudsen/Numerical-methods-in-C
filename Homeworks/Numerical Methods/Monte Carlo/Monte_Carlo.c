#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Monte_Carlo.h"
#include "Integrator.h"
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

double adapt_2D(double old_f(double x, double y),double c_f(double),
 double d_f(double), double a_old,double b_old,double acc,double eps){
	// Rescaling the points x2 and x3 and finding the function values
	int nrec=0;
	double f2, f3,a,b, result;

	double f_outer(double y) {
				double f_inner(double x){
					return old_f(x,y);
				}
				nrec=0;
				a=c_f(y); 
				b=d_f(y);
				f2=f_inner(a+2*(b-a)/6);
				f3=f_inner(a+4*(b-a)/6);
				result=adapt24(f_inner,old_f,a,b,acc,eps,f2,f3,nrec);
				return result;
			}

	a=a_old; b=b_old;
	f2=f_outer(a+2*(b-a)/6);

	f3=f_outer(a+4*(b-a)/6);
	result = adapt24(f_outer,old_f,a,b,acc,eps,f2,f3,nrec);
	return result;

}
	

