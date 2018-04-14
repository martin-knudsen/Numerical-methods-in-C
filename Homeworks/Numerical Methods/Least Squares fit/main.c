#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>9
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

double funs2(int i, double x){
   switch(i){
   case 0: return -exp(-x); break;
   case 1: return sqrt(x);   break;
   case 2: return sin(x);     break;
   case 3: return 1.0;   break;
   default: {fprintf(stderr,"funs: wrong i:%d",i); return NAN;}
   }
}

double fit(int m, gsl_vector* c, double x) {
	double sum = 0;
	for(int i=0;i<m;i++){
		sum+=gsl_vector_get(c, i)*funs2(i, x);
	}
	return sum; 
}

double fit_plus(int i, int m, gsl_vector* c, gsl_vector* dc, double x){
	double result = 0;
	result = fit(m, c, x)+gsl_vector_get(dc, i)*funs2(i,x); 
	return result;
}

double fit_minus(int i, int m, gsl_vector* c, gsl_vector* dc, double x){
	double result = 0;
	result = fit(m, c, x)-gsl_vector_get(dc, i)*funs2(i,x); 
	return result;
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

	// B.2
	
	m = 4;
	gsl_vector* c2 = gsl_vector_alloc(m);
	gsl_vector* dc2 = gsl_vector_alloc(m);
	gsl_matrix* COV2 = gsl_matrix_alloc(m, m);
	least_squares(x,y,dy,m,funs2,c2,COV2);

	for(int i=0; i<m; i++){
		gsl_vector_set(dc2, i, sqrt(gsl_matrix_get(COV2, i, i)));
	}

	
	FILE* data =fopen("data.txt", "w+");
	for(int i=0; i<n; i++){
		fprintf(data, "%g\t %g\t %g \n",ax[i], ay[i], ady[i]);
	}

	fclose(data);  

	// i=0
	FILE* fit_c0 =fopen("fit_c0.txt", "w+");
	
	double fity, fitplus, fitminus; 
	double zmax = 10, deltaz = zmax/200, z;
	
	for(z=0.0;z<zmax; z +=deltaz){
		fity = fit(m,c2,z);
		fitplus = fit_plus(0,m,c2,dc2,z);
		fitminus = fit_minus(0,m,c2,dc2,z);
		fprintf(fit_c0, "%g\t%g\t%g\t%g\n",z,fity,fitplus,fitminus);
	}
	fclose(fit_c0); 

	// i=1
	FILE* fit_c1 =fopen("fit_c1.txt", "w+");
	for(z=0.0;z<zmax; z +=deltaz){
		fity = fit(m,c2,z);
		fitplus = fit_plus(1,m,c2,dc2,z);
		fitminus = fit_minus(1,m,c2,dc2,z);
		fprintf(fit_c1, "%g\t%g\t%g\t%g\n",z,fity,fitplus,fitminus);
	}
	fclose(fit_c1); 

	// i=2
	FILE* fit_c2 =fopen("fit_c2.txt", "w+");
	for(z=0.0;z<zmax; z +=deltaz){
		fity = fit(m,c2,z);
		fitplus = fit_plus(2,m,c2,dc2,z);
		fitminus = fit_minus(2,m,c2,dc2,z);
		fprintf(fit_c2, "%g\t%g\t%g\t%g\n",z,fity,fitplus,fitminus);
	}
	fclose(fit_c2); 

	// i=3
	FILE* fit_c3 =fopen("fit_c3.txt", "w+");
	for(z=0.0;z<zmax; z +=deltaz){
		fity = fit(m,c2,z);
		fitplus = fit_plus(3,m,c2,dc2,z);
		fitminus = fit_minus(3,m,c2,dc2,z);
		fprintf(fit_c3, "%g\t%g\t%g\t%g\n",z,fity,fitplus,fitminus);
	}
	fclose(fit_c3); 

	printf("As one can see the least squares fit is not perfect but reaonably within the errorbars\n");
	printf("at least for this example purpose..\n");


	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_matrix_free(COV);
	gsl_vector_free(c2);
	gsl_vector_free(dc2);
	gsl_matrix_free(COV2);

	return EXIT_SUCCESS;
}