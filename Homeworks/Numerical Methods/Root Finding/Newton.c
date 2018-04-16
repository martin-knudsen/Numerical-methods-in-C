#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "QR.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f"

void newton(
	void f(gsl_vector* x, gsl_vector* fx),
	gsl_vector* xstart,
	double dx,
	double epsilon
){
	
}