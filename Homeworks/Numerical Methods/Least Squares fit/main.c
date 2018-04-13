#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "QR.h"
#include "QR_ls.h"
#define RND (double)rand()/RAND_MAX
#define FMT "%7.3f" //format of print "7 width, 3 digits after comma" 


// inspired by Dmitri Fedorovs print implementation
void printm(gsl_matrix *A){
	// iterating over all rows and columns 
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++) {
			printf(FMT,gsl_matrix_get(A,i,j));
		}
		printf("\n");}
}

void printv(gsl_vector *A){
	for(int i=0;i<A->size;i++){
		printf(FMT,gsl_vector_get(A,i));
		printf("\n");
	}
}

double funs(int i, double x){
   switch(i){
   case 0: return log(x); break;
   case 1: return 1.0;   break;
   case 2: return x;     break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

int main() {
		
	double ax[]  =  {0.1 ,   1.33  ,  2.55  ,  3.78   ,    5.0  ,  6.22  ,  7.45  ,  8.68  ,   9.9};
	double ay[]  =   {-15.3 ,   0.32  ,  2.45  ,  2.75  ,  2.27  ,  1.35  , 0.157  , -1.23  , -2.75};
	double ady[] =   { 1.04  , 0.594 ,  0.983  , 0.998 ,   1.11  , 0.398 ,  0.535  , 0.968  , 0.478};

	int n = sizeof(ax)/sizeof(ax[0]);
	int m = 3;
	gsl_vector* x = gsl_vector_alloc(n);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* dy = gsl_vector_alloc(n);
	gsl_vector* c = gsl_vector_alloc(m);
	gsl_matrix* COV = gsl_matrix_alloc(m, m);



	for(int i=0;i<n;i++){
		gsl_vector_set(x,i,ax[i]);
		gsl_vector_set(y,i,ay[i]);
		gsl_vector_set(dy,i,ady[i]);
	}

	least_squares(x,y,dy,m,funs,c,COV);

	printf("1.2 The value of the fit coefficients c:\n");
	printv(c);


	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_matrix_free(COV);

	return EXIT_SUCCESS;
}