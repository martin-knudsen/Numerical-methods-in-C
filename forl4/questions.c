#include <stdio.h>

int main(void)
{
	double x=1.23;
	double* p=&x;
	double x2 = *p;

	printf("%g\n",x2);

	
	
	int i=1; f1(i); printf("i=%i\n",i);
	

	
	int ii=1; f2(&ii); printf("i=%i\n",ii);
	

	
	int iii=1; f3(&iii); printf("i=%i\n",iii);
	
	return 0; 
	
}

void f1(int i){i=0;}

void f2(int* i){*i=0;}

void f3(int* i){i=NULL;}
