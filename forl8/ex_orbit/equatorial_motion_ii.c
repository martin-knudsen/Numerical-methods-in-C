#include <stdio.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_errno.h>

int 
orbital_equation (double t, const double y[], double dydt[], void * params) {
	
	(void)(t); 
	
	dydt[0] = y[1];
	dydt[1] = 1-y[0];

	return GSL_SUCCESS;
}

int main() {

	
	gsl_odeiv2_system sys;
	sys.function = orbital_equation;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = NULL;

	double hstart = 1e-3; double epsabs = 1e-6; double epsrel = 1e-6; 

	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);

	double t = 0.0; double phi_max = 20*M_PI; double phi_delta = 0.01; 
	double y[2]= {1.0, 0.5};

	/*
	printf("-----------------------------------------------\n");
	printf("Now solving the equatorial equation for eps = 0\n");
	printf("Resulting in the differential equation u''=1-u\n");
	printf("Using the starting conditions: u(0)=1 and u'(0)=0.5");
	 */
	
	for (double phi = 0; phi<=phi_max; phi += phi_delta) {
		
		int status = gsl_odeiv2_driver_apply (d, &t, phi, y);
		if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%i\n", status);
          break;
        }

      	printf ("%.5e %.5e \n", phi, y[0]);
	}
	gsl_odeiv2_driver_free (d);

	
	return 0;
}