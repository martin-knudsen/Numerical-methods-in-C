#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

typedef struct{int n; double *x,*y,*b,*c,*d;} cspline;

int search(cspline* s, double z) {
	int i=0,j=s->n-1;
	//binarysearch:
	while(j-i>1){int m=(i+j)/2;if(z>s->x[m]) i=m;else j=m;}
	return i;
}

/* Taken from Dmitri Fedorovs Interpolation chapter */
cspline* cspline_alloc(int n, double *x,double *y){//buildsnaturalcspline
	cspline* s=(cspline*) malloc(sizeof(cspline));
	s->x=(double*) malloc(n*sizeof(double));
	s->y=(double*) malloc(n*sizeof(double));
	s->b=(double*) malloc(n*sizeof(double));
	s->c=(double*) malloc((n-1)*sizeof(double));
	s->d=(double*) malloc((n-1)*sizeof(double));
	s->n=n;for(int i=0;i<n;i++){s->x[i]=x[i];s->y[i]=y[i];}
	double h[n-1],p[n-1];//VLA
	for(int i=0;i<n-1;i++){h[i]=x[i+1]-x[i];assert(h[i]>0);}
	for(int i=0;i<n-1;i++)p[i]=(y[i+1]-y[i])/h[i];
	double D[n],Q[n-1],B[n];//buildingthetridiagonalsystem:
	D[0]=2;for(int i=0;i<n-2;i++)D[i+1]=2*h[i]/h[i+1]+2;D[n-1]=2;
	Q[0]=1;for(int i=0;i<n-2;i++)Q[i+1]=h[i]/h[i+1];
	for(int i=0;i<n-2;i++)B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
	B[0]=3*p[0];B[n-1]=3*p[n-2];//Gausselimination:
	for(int i=1;i<n;i++){D[i]-=Q[i-1]/D[i-1];B[i]-=B[i-1]/D[i-1];}
	s->b[n-1]=B[n-1]/D[n-1];//back-substitution:
	for(int i=n-2;i>=0;i--)s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];
	for(int i=0;i<n-1;i++){
	s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
	s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
}
return s;
}
double cspline_eval(cspline *s,double z){
	assert(z>=s->x[0]&&z<=s->x[s->n-1]);

	int i = search(s, z);
	int j = i+1;

	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}
void cspline_free(cspline *s){
	free(s->x);free(s->y);free(s->b);free(s->c);free(s->d);free(s);}


double cspline_derivative(cspline *s, double z){
	assert(z>=s->x[0]&&z<=s->x[s->n-1]);

	int i = search(s, z);
	int j = i+1;

	double result = s->b[i]+2*s->c[i]*(z-s->x[i] \
		+s->d[i])*(3*z*z-6*z*s->x[i]+3*s->x[i]*s->x[i]);
	return result;
}
double cspline_integral(cspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);

	double integral = 0.0;

	int i_z = search(s, z);
	int j_z = i_z+1;
	
	double integ_high, integ_low, integ_ij;
	for(int i=0; i<i_z; i++) {
		int j=i+1;
		
		integ_high = s->y[i]*s->x[j]+s->b[i]*s->x[j]*(s->x[j]/2-s->x[i])-pow(s->c[i]*(s->x[i]-s->x[j]),3)/3+pow(s->d[i]*(s->x[j]-s->x[i]),4)/4;
		integ_low = s->y[i]*s->x[i]+s->b[i]*s->x[i]*(s->x[i]/2-s->x[i])-pow(s->c[i]*(s->x[i]-s->x[i]),3)/3+pow(s->d[i]*(s->x[i]-s->x[i]),4)/4;
		integ_ij = integ_high - integ_low;
		integral += integ_ij;
	}

	integ_high = s->y[i_z]*z+s->b[i_z]*z*(z/2-s->x[i_z])-pow(s->c[i_z]*(s->x[i_z]-z),3)/3+pow(s->d[i_z]*(z-s->x[i_z]),4)/4;
	integ_low = s->y[i_z]*s->x[i_z]+s->b[i_z]*s->x[i_z]*(s->x[i_z]/2-s->x[i_z])-pow(s->c[i_z]*(s->x[i_z]-s->x[i_z]),3)/3+pow(s->d[i_z]*(s->x[i_z]-s->x[i_z]),4)/4;
	integ_ij = integ_high - integ_low;
	integral += integ_ij;

	return integral;
}

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

	FILE* output = fopen("output.txt", "a");

	fprintf(output,"Testing if cspline spline works. \n");
	fprintf(output,"Testfunction is cos(x) from 0 to 2pi. \n");
	fprintf(output,"Data for plotting the cspline is found in cdata.txt \n");
	
	FILE *cdata = fopen("cdata.txt", "w+");

	double delta_x = 2*M_PI/100, x_start = 0.0, x_slut=2*M_PI, x_it;
	double spline_result;

	cspline* cspline_cos = cspline_alloc(n, x, y);
	for(x_it = x_start; x_it <x_slut; x_it += delta_x) {
		
		spline_result = cspline_eval(cspline_cos, x_it);
		fprintf(cdata, "%g \t %g \n", x_it, spline_result);
	}

	fclose(cdata);
	
	fprintf(output,"Integral of cos(x) from 0 to pi/2 is 1 analytically \n");

	double integ_c  = cspline_integral(cspline_cos, M_PI/2);
	fprintf(output,"Using cubic splines I get %g\n", integ_c);

	fprintf(output,"The theoretical derivative at point pi/2 is -1\n");

	double der_c  = cspline_derivative(cspline_cos, M_PI/2);
	fprintf(output,"Using cubic splines I get %g \n", der_c); 
	

	cspline_free(cspline_cos);

	/* extra test using gsl spline and interp*/
	const gsl_interp_type* T = gsl_interp_cspline;
	gsl_spline* cspline_gsl = gsl_spline_alloc(T, n);

	int status = gsl_spline_init(cspline_gsl, x, y, n);

	gsl_interp_accel *acc = gsl_interp_accel_alloc();

	FILE* csplinegsldata = fopen("csplinegsldata.txt","w+");
	for(double z=0.0; z<=2*M_PI;z += 2*M_PI/100) {

		spline_result = gsl_spline_eval(cspline_gsl, z, acc);
		fprintf(csplinegsldata, "%g \t %g \n", z, spline_result);
	}

	fclose(csplinegsldata);
	fprintf(output,"I have plottet the result from gsl and my cspline in csplinetest.pdf\n");
	void gsl_spline_free(cspline_gsl);

	double derivgsl = gsl_spline_eval_deriv(cspline_gsl, M_PI/2, acc);
	fprintf(output,"This is the derivative result at pi/2 from GSL: %g\n", derivgsl);
	double integgsl = gsl_spline_eval_integ(cspline_gsl, 0, M_PI/2, acc);
	fprintf(output,"This is the integral result at pi/2 from GSL: %g\n", integgsl);
	fprintf(output,"\n");
	fclose(output);
	return EXIT_SUCCESS;
}
