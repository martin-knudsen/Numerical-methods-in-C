#include <stdio.h>
#include <math.h>
#include <gsl/gsl_odeiv2.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

int 
error_function(double t, const double y[], \
	double dydt[], void * params) {

	/* because they are not used explicitly */
	(void)(params);
	(void)(y);

	dydt[0]=2/sqrt(M_PI)*exp(-t*t);

	return GSL_SUCCESS;
}


void master(double x_start, double x_end, double dx) {

	gsl_odeiv2_system system;

	system.function = error_function;
	system.dimension = 1;
	system.jacobian = NULL;
	system.params = NULL;

	

	double hstart = 1e-3, epsabs = 1e-6, epsrel = 1e-6;

	gsl_odeiv2_driver *driver = \
	gsl_odeiv2_driver_alloc_y_new(\
		&system, \
		gsl_odeiv2_step_rk8pd, \
		hstart,\ 
		epsabs,\
		epsrel); 

	double x = x_start, y[1]={0.0};

	for (double xx = 0.0; xx<=x_end; xx += dx) {
		
		int status = gsl_odeiv2_driver_apply (driver, &x, xx, y);
		if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%i\n", status);
          break;
        }

      	printf ("%.5e \t %.5e \n", xx, y[0]);
	}
	gsl_odeiv2_driver_free (driver);

}

int main(int argc, const char** argv) {
	
	double x_start = atof(argv[1]);
	double x_end = atof(argv[2]);
	double dx = atof(argv[3]);

	master(x_start, x_end, dx);
	return EXIT_SUCCESS;
}