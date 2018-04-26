#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

/* These methods are taken from Dmitri Fedorovs lecture notes 
 It works by taking in the function, integration limits a and b, the accuracies, 
 and the previous calculated values of the function at x2 and x3. It outputs
 the value of the integral evaluated on the limits within the tolerances, and if 
 it is not within toleranced it recures until it has found the answer. 
 */
double adapt24(double f(double),double a,double b,
	double acc,double eps,double f2,double f3,int nrec){
	// upper limit on how many recursions are possible until if will just accept the error
	assert(nrec<1000000);
	//Rescaling to the original interval point x1 and x4 and finding the function value
	double f1=f(a+(b-a)/6),f4=f(a+5*(b-a)/6);
	// calculate integral according to Q, the trapez method and q the rectangular method
	// according to (8.49) and (8.50). The weights have also been rescaled by (b-a).
	double Q=(2*f1+f2+f3+2*f4)/6*(b-a),q=(f1+f4+f2+f3)/4*(b-a);
	// calculate the tolerance and error according to (8.46), (8.47)
	double tolerance=acc+eps*fabs(Q),error=fabs(Q-q);
	// See if the error is acceptable then can return the estimation of the integral
	// of the subinterval
	if(error<tolerance)return Q;
	// if it isn't acceptable, subdivide the interval into to equally big intervals
	// update the absolute tolerance and try again on each new subinterval. This is the recursion
	// part
	else{
		double Q1=adapt24(f,a,(a+b)/2,acc/sqrt(2.),eps,f1,f2,nrec+1);
		double Q2=adapt24(f,(a+b)/2,b,acc/sqrt(2.),eps,f3,f4,nrec+1);
	// return the sum of the integral of both subintervals
	return Q1+Q2;}
}

// The first step of calculating the function values at f2 and f3
double adapt(
	double f(double),double a,double b,double acc,double eps){
	// Rescaling the points x2 and x3 and finding the function values
	double f2=f(a+2*(b-a)/6),f3=f(a+4*(b-a)/6);int nrec=0;
	return adapt24(f,a,b,acc,eps,f2,f3,nrec);
}
