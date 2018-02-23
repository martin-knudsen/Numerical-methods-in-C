#include <stdio.h>
#include <gsl/gsl_sf_airy.h>

int () {

	double x;
	while(scanf("%lg",&x) != EOF) {
		double Ai = gsl_sf_airy_Ai(x);
		double Bi = gsl_sf_airy_Bi(x);
		printf("%lg \t %g \t %g\n", x, Ai, Bi);
	}
	return 0;
}
