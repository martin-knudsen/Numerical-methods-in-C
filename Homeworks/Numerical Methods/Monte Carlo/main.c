#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include "Integrator.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.6f" //format of print "7 width, 3 digits after comma" 

int main(void) {
	// part A 
	printf("A. Plain Monte Carlo integration\n\n");
	printf("Integrating dr dphi dtheta r^2 from 0 to 1, 0 to 2pi and 0 to pi\n");
	
	int dim=3, N=1e6; double a[3]={0,0,0}, b[3]={1,2*M_PI,M_PI},result, error;
	printf("Number of points: %i\n",N);
	double f1(double* x){
		return x[0]*x[0];
		
	}

	plainmc(dim,a,b,f1,N,&result,&error);

	printf("Result: %g, correct answer is: %g\n",result,2*M_PI*M_PI/3);
	printf("Error: %g \n\n",error);

	printf("Integrating dx dy dz dq x*y^2*z^3*q^4 from 0 to 1, 1 to 2 and 2 to 3 and 3 to 4\n");
	
	dim=4; N=1e7; double a2[4]={0,1,2,3}, b2[4]={1,2,3,4};
	printf("Number of points: %i\n",N);
	double f2(double* x){
		return x[0]*x[1]*x[1]*x[2]*x[2]*x[2]*x[3]*x[3]*x[3]*x[3];
	}

	plainmc(dim,a2,b2,f2,N,&result,&error);

	printf("Result: %g, correct answer is: %g\n",result,71071.0/24.0);
	printf("Error: %g\n\n",error);

	printf("Integrating dx dy dz (1-cos(x)cos(y)cos(z))^-1/pi^3 from 0 to pi, 0 to pi and 0 to pi\n");
	
	dim=3; N=1e7; double a3[3]={0,0,0}, b3[3]={M_PI,M_PI,M_PI};
	printf("Number of points: %i\n",N);
	double f3(double* x){
		return 1/(1-cos(x[0])*cos(x[1])*cos(x[2]))/(M_PI*M_PI*M_PI);
	}

	plainmc(dim,a3,b3,f3,N,&result,&error);

	printf("Result: %g, correct answer is: %7.6f \n",result,1.3932039296856768591842462603255);
	printf("Error: %g\n\n",error);


	// part B
	printf("B. Check that error behaves as O(1/√N)\n\n");	

	// part C
	printf("C. 2D adaptive integrator\n\n");


	return EXIT_SUCCESS;
}