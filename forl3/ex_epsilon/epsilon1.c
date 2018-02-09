#include <limits.h>
#include <float.h>
#include <stdio.h>
#include "tau.h"
int main()
{
	printf("Exercise 1.i\n");

	printf("The real int_max is: %i \n",INT_MAX);
	printf("The while loop int_max is: %i\n",int_max_while());
	printf("The for loop int_max is: %i\n",int_max_for());
	printf("The do while loop int_max is: %i\n",int_max_do_while());

	printf("Exercise 1.ii\n");

	printf("The real int_min is: %i \n",INT_MIN);
	printf("The while loop int_min is: %i\n",int_min_while());
	printf("The for loop int_min is: %i\n",int_min_for());
	printf("The do while loop int_min is: %i\n",int_min_do_while());

	printf("Exercise 1.iii\n");
	
	printf("Real epsilon float is: %f \n",FLT_EPSILON);

	float f1=1.0; 
	while(1+f1!=1)
		{
			f1/=2;
		} 
	f1*=2;
	printf("While loop is: %f \n",f1);
		
	

	float f2=0;
	for(float e=1; 1+e!=1; e/=2)
	{f2=e;}
	printf("For loop is: %f \n", f2);
		
	float f3 =1;
	do  
	{
		f3/=2;

	}
	while(f3+1!=1.0);
	f3*=2;
	printf("Do while loop is: %f \n", f3);

	printf("Real epsilon double is: %g \n",DBL_EPSILON);
	double d1=1; 
	while(1+d1!=1)
		{
			d1/=2;
		} 
	d1*=2;
	printf("While loop is: %g \n",d1);

	double d2 =0;
	for(double e=1; 1+e!=1; e/=2)
	{d2=e;}
	printf("For loop is: %g \n",d2);
		

	double d3 =1;
	do  
	{
		d3/=2;

	}
	while(d3+1!=1.0);
	d3*=2;
	printf("Do while loop is: %g \n",d3);
	printf("Real epsilon long double is: %LG \n",LDBL_EPSILON);

	long double ld1=1; 
	while(1+ld1!=1)
		{
			ld1/=2;
		} 
	ld1*=2;
	printf("While loop is: %LG \n",ld1);
	
	long double ld2=0;
	for(long double e=1; 1+e!=1; e/=2)
	{ld2 = e;}
	printf("For loop is: %LG \n",ld2);
	
	long double ld3 =1;
	do  
	{
		ld3/=2;

	}
	while(1+ld3!=1.0);
	ld3*=2;
	printf("Do while loop is: %LG \n",ld3);


	printf("Exercise 2.i\n");

	int max=INT_MAX/2;

	float sum_up_float = 0.0;
	for(int i =1;i<max+1;i++)
	{
		sum_up_float += 1.0f/i;
	}

	printf("sum_up_float has value=%f\n", sum_up_float);

	float sum_down_float = 0.0;
	for(int i =max;i>0;i--)
	{
		sum_down_float += 1.0f/i;
	}
	printf("sum_down_float has value=%f\n", sum_down_float);
	printf("Exercise 2.ii\n");
	printf("The two functions are mirrors of eachother, so the only difference is\
		that one starts at max and iterates down and the other opposite.\
		however we observe a difference which must be attributed to this not being precise\n");

	printf("Exercise 2.iii\n");

	printf("Yes it converges towards 1.0 because\
		then all the terms excluding 1.0 are zero\n");

	printf("Exercise 2.iv\n");

	double sum_up_double = 0.0;
	for(int i =1;i<max+1;i++)
	{
		sum_up_double += 1.0/i;
	}

	printf("sum_up_double has value=%g\n", sum_up_double);

	double sum_down_double = 0.0;
	for(int i =max;i>0;i--)
	{
		sum_down_double += 1.0/i;
	}
	printf("sum_down_double has value=%g\n", sum_down_double);

	printf("We observe a clear difference in the summing of the float\
		and double. Somehow the float is much more sensitive than\
		the double.");
	

	return 0;
}

int int_max_while()
{
	int i=0;
	while(i+1>i)
	{
		i++;
	}
	return i;
}

int int_max_for()
{
	int int_max=0;
	for(int i;i+1<i;i++)
	{
		int_max = i;
	}
	return int_max;
}

int int_max_do_while()
{
	int i=0;
	do 
	{
		i++;
	}
	while(i+1>i);
	return i;
}

int int_min_while()
{
	int i=0;
	while(i-1<i)
	{
		i--;
	}
	return i;
}

int int_min_for()
{
	int int_min=0;
	for(int i;i-1>=i;i--)
	{
		int_min = i;
	}
	return int_min;
}

int int_min_do_while()
{
	int i=0;
	do 
	{
		i--;
	}
	while(i-1<i);
	return i;
}