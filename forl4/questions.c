#include <stdio.h>
int i=2;

void f(){printf("i=%i\n",i);}
int main(void)
{
	double x=1.23;
	double* p=&x;
	double x2 = *p;

	printf("%p\n",p);

	
	
	int i2=1; f1(i2); printf("i=%i\n",i2);
	

	
	int ii=1; f2(&ii); printf("i=%i\n",ii);
	

	
	int iii=1; f3(&iii); printf("i=%i\n",iii);




	int i=1; /* function scope */
	{
		int i=0; /* block scope */
		printf("i=%i\n",i);
	}
	printf("i=%i\n",i);
	f();
	
	
	return 0; 
	
}

void f1(int i2){i2=0;}

void f2(int* i2){*i2=0;}

void f3(int* i2){i2=NULL;}
