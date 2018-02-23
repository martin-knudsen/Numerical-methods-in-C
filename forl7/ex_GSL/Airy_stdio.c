#include <stdio.h>
#include <gsl/gsl_sf_airy.h>

int main() {

	double x;
	while(scanf("%lg",&x) != EOF) {
		double Ai = gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
		double Bi = gsl_sf_airy_Bi(x, GSL_PREC_DOUBLE);
		printf("%lg \t %g \t %g\n", x, Ai, Bi);
	}
	return 0;
}
