#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <math.h>

double
test_function (double x, void * p) {
   (void)(p);
   
   return  log(x)/sqrt(x);
}



int main() {
	
	gsl_function F;
	F.function = test_function;
	F.params = NULL;
	
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);

    double expected_result = -4.0;
    double a = 0.0; double b= 1.0; 
	double epsabs = 1e-6; double epsrel = 1e-6; double limit= 1000;
	double result, abserr;

  gsl_integration_qags (&F, a, b, epsabs, epsrel, limit, workspace, &result, &abserr);

  printf("Integrating ln(x)/sqrt(x) from 0 to 1\n");
  printf("Cecause it contains a singularity at zero we use qaqs\n");
  printf ("Result = % g\n", result);
  printf ("Expected result = % g\n", expected_result);

  gsl_integration_workspace_free (workspace);
	
	return 0; 
}