#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <math.h>
#include "normintegral.h"
#include "hamiltonintegral.h"

int main() {

	double alfa_start = 0.01;
	double alfa_max = 4.00;
	double alfa_delta = 0.01;


	for(double alfa = alfa_start; alfa<alfa_max; alfa += alfa_delta) {
		double E_alfa = integrate_ham(&alfa)/integrate_norm(&alfa);
		printf("%g \t %g\n", alfa, E_alfa);
	}
	
	
	return 0; 
}