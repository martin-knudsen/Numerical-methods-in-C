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