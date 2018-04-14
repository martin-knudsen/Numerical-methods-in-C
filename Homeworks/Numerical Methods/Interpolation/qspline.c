#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

/* Taken from Dmitri Fedorovs Interpolation chapter,
   explained in the linearspline.c program.*/
int search(qspline* s, double z) {
	int i=0,j=s->n-1;
	//binarysearch:
	while(j-i>1){int m=(i+j)/2;if(z>s->x[m]) i=m;else j=m;}
	return i;
}

/*Also taken from Dmitri. We use the special struct saving x and y lists
as well as all the coefficients needed and the number of data points.
 */
qspline* qspline_alloc(int n, double* x,double* y){//buildsqspline
	// start by allocating everything for the struct
	qspline *s=(qspline*) malloc(sizeof(qspline));//spline
	s->b=(double*) malloc((n-1)*sizeof(double)); //bi
	s->c=(double*) malloc((n-1)*sizeof(double)); //ci
	s->x=(double*) malloc(n*sizeof(double)); //xi
	s->y=(double*) malloc(n*sizeof(double)); //yi
	s->n=n; 
	// transfer x and y to struct
	for(int i=0; i<n; i++){s->x[i]=x[i];s->y[i]=y[i];}
	
	// initiate and make helper list h and coefficient list p according
	// to (1.6)
	int i; double p[n-1], h[n-1]; //VLAfromC99
	for(i=0; i<n-1; i++){h[i]=x[i+1]-x[i];p[i]=(y[i+1]-y[i])/h[i];}
	
	// perform first recursion for coefficient c upwards and directly into struct
	s->c[0]=0; //recursionup:
	for(i=0;i<n-2;i++) s->c[i+1]=(p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];
	
	// perform second for downwards with proper starting point into struct
	s->c[n-2]/=2; //recursiondown:
	for(i=n-3;i>=0;i--) s->c[i]=(p[i+1]-p[i]-s->c[i+1]*h[i+1])/h[i];
	// these two followed (1.11) and (1.12) in lecture notes

	// calculate b's from c and p because more convenient
	// add to struct
	for(i=0;i<n-1;i++) s->b[i]=p[i]-s->c[i]*h[i];
	return s;}

double qspline_eval(qspline *s, double z){
	//evaluatess(z)
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = search(s, z);
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);
}
//interpolatingpolynomial

void qspline_free(qspline *s){//freetheallocatedmemory
	free(s->x); 
	free(s->y); 
	free(s->b);
	free(s->c);
	free(s);
}

double qspline_derivative(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	
	int i = search(s, z);
	
	double result = s->b[i]+2*s->c[i]*(z-s->x[i]);
	return result;
} /* evaluates the derivative of the prebuilt spline at point z */

double qspline_integral(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);

	double integral = 0;

	int i_z = search(s, z);
	int j_z = i_z+1;
	double integ_high, integ_low, integ_ij;
	for(int i=0; i<i_z; i++) {
		int j=i+1;

		integ_high = s->y[i]*s->x[j]+s->b[i]*s->x[j]*(0.5*s->x[j]-s->x[i])+\
		s->c[i]*s->x[j]*(1.0/3.0*s->x[j]*s->x[j]+s->x[i]*s->x[i]-s->x[j]*s->x[i]);
		integ_low = s->y[i]*s->x[i]+s->b[i]*s->x[i]*(0.5*s->x[i]-s->x[i])+\
		s->c[i]*s->x[i]*(1.0/3.0*s->x[i]*s->x[i]+s->x[i]*s->x[i]-s->x[i]*s->x[i]); 
		integ_ij = integ_high - integ_low;
		integral += integ_ij;
	}
	
	integ_high = s->y[i_z]*z+s->b[i_z]*z*(0.5*z-s->x[i_z])+\
	s->c[i_z]*z*(1.0/3.0*z*z+s->x[i_z]*s->x[i_z]-z*s->x[i_z]);
	integ_low = s->y[i_z]*s->x[i_z]+s->b[i_z]*s->x[i_z]*(0.5*s->x[i_z]-s->x[i_z])+\
	s->c[i_z]*s->x[i_z]*(1.0/3.0*s->x[i_z]*s->x[i_z]+s->x[i_z]*s->x[i_z]-s->x[i_z]*s->x[i_z]); 
	integ_ij = integ_high - integ_low;
	integral += integ_ij;
	return integral;
} /* evaluates the integral of the prebuilt spline from x[0] to z */

/* function used is cosine from 0 to 2pi using 11 points */
int main() {
	
	const int n = 11;  
	double x[11] = {0.0, 2*M_PI/20*2, 2*M_PI/20*4, \
		2*M_PI/20*6, 2*M_PI/20*8, \
		2*M_PI/20*10, 2*M_PI/20*12,  \
		2*M_PI/20*14, 2*M_PI/20*16, \
		2*M_PI/20*18, 2*M_PI/20*20};
	double y[11] = {1.0, 0.80901,  0.30901, \
		 -0.30901,   -0.80901,   -1.0, \
		   -0.80901,    -0.30901,  \
		 0.30901, 0.80901,  1.0}; 

	FILE* output = fopen("out.txt", "a");

	fprintf(output,"Testing if qspline spline works. \n");
	fprintf(output,"Testfunction is cos(x) from 0 to 2pi. \n");
	fprintf(output,"Data for plotting the qspline is found in qdata.txt \n");
	
	FILE *qdata = fopen("qdata.txt", "w+");

	double delta_x = 2*M_PI/100, x_start = 0.0, x_slut=2*M_PI, x_it;
	double spline_result;

	qspline* qspline_cos = qspline_alloc(n, x, y);
	for(x_it = x_start; x_it <x_slut; x_it += delta_x) {
		
		spline_result = qspline_eval(qspline_cos, x_it);
		fprintf(qdata, "%g \t %g \n", x_it, spline_result);
	}

	fclose(qdata);
	
	fprintf(output,"Integral of cos(x) from 0 to pi/2 is 1 analytically \n");

	double integ_q  = qspline_integral(qspline_cos, M_PI/2);
	fprintf(output,"Using qaudric splines I get %g\n", integ_q);

	fprintf(output,"The theoretical derivative at point pi/2 is -1\n");

	double der_q  = qspline_derivative(qspline_cos, M_PI/2);
	fprintf(output,"Using quadric splines I get %g \n", der_q); 
	

	qspline_free(qspline_cos);

	/* extra test */
	int n_test = 5;
	// data
	double x1[5] = {1, 2, 3, 4, 5};
	double x2[5] = {1, 2, 3, 4, 5};
	double x3[5] = {1, 2, 3, 4, 5};

	double y1[5] = {1, 1, 1, 1, 1};
	double y2[5] = {1, 2, 3, 4, 5};
	double y3[5] = {1, 4, 9, 16, 25};

	//predicted values of coefficients
	double p1[4] = {0, 0, 0, 0};
	double p2[4] = {1, 1, 1, 1};
	double p3[4] = {3, 5, 7, 9};

	double c1[4] = {0, 0, 0, 0};
	double c2[4] = {0, 0, 0, 0};
	double c3[4] = {1, 1, 1, 1};

	double b1[4] = {0, 0, 0, 0};
	double b2[4] = {1, 1, 1, 1};
	double b3[4] = {2, 4, 6, 8};
	
	// spline objects

	qspline* test1 = qspline_alloc(n_test, x1, y1);
	qspline* test2 = qspline_alloc(n_test, x2, y2);
	qspline* test3 = qspline_alloc(n_test, x3, y3);

	fprintf(output,"Testing that all the coefficients calculated manually are the same\n");
	fprintf(output,"As I get using the spline program\n");

	bool status = true;

	for(int i=0; i<n_test-1; i++) {
		if(c1[i]!=test1->c[i]) status = false;
		if(c2[i]!=test2->c[i]) status = false;
		if(c3[i]!=test3->c[i]) status = false;
		if(b1[i]!=test1->b[i]) status = false;
		if(b2[i]!=test2->b[i]) status = false;
		if(b3[i]!=test3->b[i]) status = false;
	} 
	fprintf(output,"The result of the test is: \n");
	fprintf(output,status ? "True\n" : "False\n");
	fprintf(output,"\n");
	fclose(output);
	return EXIT_SUCCESS;
}
