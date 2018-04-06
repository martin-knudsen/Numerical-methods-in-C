#include <stdio.h>
#include "komplex.h"

int main(){
	komplex a = {1,2}, b = {3,4}, c = {1,2};
	komplex  z;

	printf("testing komplex_set...\n");
	komplex_print("a=", a);
	printf("setting z as a\n");
	komplex_set(&z, a.re, a.im);
	printf("Now printing z with komplex_print\n");
	komplex_print("z=", z);

	printf("Now testing komplex_add...\n");
	komplex_set(&a,a.re, a.im);
	komplex_set(&b,b.re, b.im);
	komplex_set(&c,c.re, c.im);

	komplex a_b_sum = komplex_add(a, b);
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("a+b=",a_b_sum);

	printf("Now testing komplex_equal...\n");
	int a_b_equal = komplex_equal(a,b);
	int a_c_equal = komplex_equal(a,c);
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("c=", c);
	printf("(a==b)=%i\n", a_b_equal);
	printf("(a==c)=%i\n", a_c_equal);

	printf("Now testing komplex_conjugate...\n");
	komplex a_conj = komplex_conjugate(a);
	komplex_print("a=", a);
	komplex_print("a*=",a_conj);

	printf("Now testing komplex_equal...\n");
	komplex a_b_div = komplex_div(a,b);
	komplex a_c_div = komplex_div(a,c);
	komplex_print("a=", a);
	komplex_print("b=", b);
	komplex_print("c=", c);
	komplex_print("(a/b)=", a_b_div);
	komplex_print("(a/c)=", a_c_div);

	printf("Now testing komplex_abs...\n");
	komplex a_abs= komplex_abs(a);
	komplex_print("a=", a);
	komplex_print("|a|=",a_abs);

	printf("Now testing komplex_exp...\n");
	komplex a_exp= komplex_exp(a);
	komplex_print("a=", a);
	komplex_print("exp(a)=",a_exp);

	printf("Now testing komplex_sin...\n");
	komplex a_sin= komplex_sin(a);
	komplex_print("a=", a);
	komplex_print("sin(a)=",a_sin);

	printf("Now testing komplex_cos...\n");
	komplex a_cos= komplex_cos(a);
	komplex_print("a=", a);
	komplex_print("cos(a)=",a_cos);

	printf("Now testing complex_sqrt...\n");
	komplex a_sqrt= komplex_sqrt(a);
	komplex_print("a=", a);
	komplex_print("sqrt(a)=+/-",a_sqrt);
	
}
