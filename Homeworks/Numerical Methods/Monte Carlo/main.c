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
	
	// part B
	printf("B. Check that error behaves as O(1/√N)\n\n");	

	// part C
	printf("B. 2D adaptive integrator\n\n");


	return EXIT_SUCCESS;
}