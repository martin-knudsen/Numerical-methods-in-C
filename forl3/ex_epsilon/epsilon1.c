#include <limits.h>
#include <float.h>
#include <stdio.h>
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
	float eps_float[3]
	float eps_float[] = epsilon_float();
	double eps_double[] = epsilon_double();
	long double eps_long_double[] = epsilon_long_double();

	printf("Real epsilon float is: %g \n",FLT_EPSILON);
	printf("While loop is: %g \n",eps_float[0]);
	printf("For loop is: %g \n",eps_float[1]);
	printf("Do while loop is: %g \n",eps_float[2]);

	printf("Real epsilon double is: %g \n",DBL_EPSILON);
	printf("While loop is: %g \n",eps_double[0]);
	printf("For loop is: %g \n",eps_double[1]);
	printf("Do while loop is: %g \n",eps_double[2]);

	printf("Real epsilon long double is: %LG \n",LDBL_EPSILON);
	printf("While loop is: %g \n",eps_long_double[0]);
	printf("For loop is: %g \n",eps_long_double[1]);
	printf("Do while loop is: %g \n",eps_long_double[2]);

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

float epsilon_float()
{
	float x=1; 
	while(1+x!=1)
		{
			x/=2;
		} 
	

	
	for(float e=1; 1+e!=1; e/=2)
	{}
	
	float a =1;
	do  
	{
		a/=2;

	}
	while(a!=1.0);

	float epsilon[] = {x, e, a}
	return epsilon
}

double epsilon_double()
{
	double x=1; 
	while(1+x!=1)
		{
			x/=2;
		} 


	for(double e=1; 1+e!=1; e/=2)
	{}
	
	double a =1;
	do  
	{
		a/=2;

	}
	while(a!=1.0);

	double epsilon[] = {x, e, a}
	return epsilon
}

long double epsilon_long_double()
{
	long double x=1; 
	while(1+x!=1)
		{
			x/=2;
		} 
	

	for(long double e=1; 1+e!=1; e/=2)
	{}
	
	long double a =1;
	do  
	{
		a/=2;

	}
	while(a!=1.0);

	long double epsilon[] = {x, e, a}
	return epsilon
}