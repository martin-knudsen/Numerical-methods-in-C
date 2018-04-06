#include <complex.h>
#include <stdio.h>
#include <math.h>
int main()
{	
	printf("Gamma (5) = %f\n, j0(0.5)=%f\n",tgamma(5), j0(0.5));
	double a_r = creal(csqrt(-2));
	double a_i = cimag(csqrt(-2));
	
	double b_r = creal(cpow(M_E,I));
	double b_i = cimag(cpow(M_E,I));

	double c_r = creal(cpow(M_E,M_PI*I));
	double c_i = cimag(cpow(M_E,M_PI*I));

	double d_r = creal(cpow(I,-M_E));
	double d_i = cimag(cpow(I,-M_E));
	
	printf("a = %f + i*%f\n", a_r, a_i);
	printf("b = %f + i*%f\n", b_r, b_i);
	printf("c = %f + i*%f\n", c_r, c_i);
	printf("d = %f + i*%f\n", d_r, d_i);


	float n1 = 0.1111111111111111111111111111;
	double n2 = 0.1111111111111111111111111111;
	long double n3 = 0.1111111111111111111111111111L;
	printf("float er %.25g\n double er %.25lg\n long double er %.25Lg",n1,n2,n3);


}

