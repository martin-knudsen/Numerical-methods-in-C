#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include "hamiltonintegral.h"


double
hamilton_function (double x, void * p) {
    double alfa  = *(double *) p;
    double f = (-alfa*alfa*x*x/2+alfa/2+x*x/2)*exp(-alfa*x*x);
  	
	return  f;
}

double integrate_ham(double* alfa) {
 	
	gsl_function F;
	F.function = hamilton_function;
	F.params = (void*)alfa;;
	
	gsl_integration_workspace * workspace = gsl_integration_workspace_alloc (1000);
 
	double epsabs = 1e-6; double epsrel = 1e-6; double limit= 1000;
	double result, abserr;

  	gsl_integration_qagi (&F, epsabs, epsrel, limit, workspace, &result, &abserr);

  	gsl_integration_workspace_free (workspace);
	
	return result; 
}