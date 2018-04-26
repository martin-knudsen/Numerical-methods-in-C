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

	double a=0,b=1,acc=1e-3,eps=1e-6; int calls=0;

	double f1(double x) {calls++; return sqrt(x);}
	double f2(double x) {calls++; return 1.0/sqrt(x);}
	double f3(double x) {calls++; return log(x)/sqrt(x);}
	double f4(double x) {calls++; return 4*sqrt(1-(1-x)*(1-x));}
	

	double result1 = adapt(f1,a,b,acc,eps); 



	printf("A. Recursive adaptive integration\n\n");
	
	printf("Integrating sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=2/3, Calculated=%g\n",result1);
	printf("Number of function calls: %i\n\n",calls); calls=0;
	double result2 = adapt(f2,a,b,acc,eps); 
	printf("Integrating 1/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=2, Calculated=%g\n",result2);
	printf("Number of function calls: %i\n\n",calls);calls=0;
	double result3 = adapt(f3,a,b,acc,eps); 
	printf("Integrating ln(x)/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=-4, Calculated=%g\n",result3);
	printf("Number of function calls: %i\n\n",calls);calls=0;
	eps=1e-18; acc=1e-18; 
	double result4 = adapt(f4,a,b,acc,eps);
	printf("Integrating 4*sqrt(1-(1-x)^2) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=pi, Calculated=%7.20f\n",result4);
	printf("Number of function calls: %i\n\n\n",calls);

	// B
	printf("B. Infinite limits\n\n");

	a=0; b=INFINITY;acc=1e-3;eps=1e-6; calls=0;

	double integral_function1(double x){
		double a = 2.0;
		int n = 1;
		calls++;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}

	double integral_function2(double x){
		double a = 2.0;
		int n = 3;
		calls++;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}

	double gaussian_theo(int n){
		double a=2.0;
		double fac=1;
		double fac2=1;
		for(int i=2;i<=n;i++){
			fac*=i;
		}
		for(int i=2;i<=2*n;i++){
			fac2*=i;
		}
		double integral = sqrt(M_PI)*fac2/fac*pow((a/2),2*n+1);
		return integral;
	}

	double integral_function1_gsl(double x, void* params){
		double a = 2.0;
		int n = 1;
		calls++;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}

	double integral_function2_gsl(double x, void* params){
		double a = 2.0;
		int n = 3;
		calls++;
		return pow(x,2*n)*exp(-x*x/(a*a));
	}


	double result5 = adapt(integral_function1,a,b,acc,eps);
	printf("Integrating x^2*exp(-x^2/4) from 0 to infinity\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=%g, Calculated=%g\n",gaussian_theo(1),result5);
	printf("Number of function calls: %i\n",calls); calls=0;


	gsl_function F;
	F.function = integral_function1_gsl;
	F.params = NULL;
	
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);
 	
	double limit= 1000;
	double abserr;
	result1=0;
  	gsl_integration_qagiu (&F, 0.0,acc, eps, limit, workspace, &result1, &abserr);
  	printf("GSL_result is:%g\n",result1);
  	printf("Number of function calls: %i\n\n",calls); calls=0;

  	gsl_integration_workspace_free (workspace);


	double result6 = adapt(integral_function2,a,b,acc,eps); 

	printf("Integrating x^6*exp(-x^2/4) from 0 to infinity\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=%g, Calculated=%g\n",gaussian_theo(3),result6);
	printf("Number of function calls: %i\n",calls);calls=0;

	result2=0;
	F.function = integral_function2_gsl;
	F.params = NULL;
	
	gsl_integration_workspace * workspace2 = gsl_integration_workspace_alloc (1000);
 	
  	gsl_integration_qagiu (&F, 0.0, acc, eps, limit, workspace2, &result1, &abserr);
  	printf("GSL_result is:%g\n",result1);
  	printf("Number of function calls: %i\n\n",calls); calls=0;

	gsl_integration_workspace_free (workspace2);

	printf("As one can see the gsl library was more effective in the second case\n"
		"but not in the first. It uses the Gauss-Kronrod 21 algorithm with \n"
		"variable transformations as I did. So perhaps their algorithm is best for complicated\n"
		"functions\n\n");

	// C
	
	printf("C. Clenshaw–Curtis variable transformation\n\n");

	a=0; b=1;acc=1e-3,eps=1e-6; calls=0;
	double result7 = clenshaw_curtis(f2,a,b,acc,eps);
	printf("Integrating 1/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=2, Calculated=%g\n",result7);
	printf("Number of function calls: %i\n",calls); calls=0;
	
	double result8 = clenshaw_curtis(f3,a,b,acc,eps);
	printf("Integrating ln(x)/sqrt(x) from 0 to 1\n");
	printf("Absolute accuracy:%g, relative accuracy: %g\n",acc,eps);
	printf("Theoretical=-4, Calculated=%g\n",result8);
	printf("Number of function calls: %i\n\n",calls); calls=0;

	printf("As one can see it is actually many times faster than the\n"
		"simple adaptive one we made before for integrations diverging at the limit.\n");

	return EXIT_SUCCESS;
}