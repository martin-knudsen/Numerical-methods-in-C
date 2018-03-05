#include"stdio.h"
#include"math.h"
double f(double x){return  A*exp(-x/T)+B;}
int main(){
	double A=3.5, T=3.2, B=1.2;
	for(double t=0;t<=5;t+=1) printf("%g %g\n",t,f(t));
return 0;
}
