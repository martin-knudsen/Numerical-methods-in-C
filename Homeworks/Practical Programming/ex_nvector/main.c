#include "nvector.h"
#include "stdio.h"
#include "stdlib.h"
#define RND (double)rand()/RAND_MAX

int main()
{
	int n = 5;

	printf("\nmain: testing nvector_alloc ...\n");
	nvector *v = nvector_alloc(n);
	if (v == NULL) printf("test failed\n");
	else printf("test passed\n");

	printf("\nmain: testing nvector_set and nvector_get ...\n");
	double value = RND;
	int i = n / 2;
	nvector_set(v, i, value);
	double vi = nvector_get(v, i);
	if (vi==value) printf("test passed\n");
	else printf("test failed\n");
	
	printf("\nmain: testing nvector_dot_product ...\n");
	nvector *a =nvector_alloc(2);
	nvector *b =nvector_alloc(2);

	nvector_set(a, 0, (double) 3.0);
	nvector_set(a, 1, (double) 4.0);
	nvector_set(b, 0, (double) 5.0);
	nvector_set(b, 1, (double) 1.0);

	double dot_product = nvector_dot_product(a,b);
	nvector_print("a =", a);
	nvector_print("b =", b);
	printf("a*b =%g\n", dot_product);

	printf("\nmain: testing nvector_set_zero ...\n");
	nvector_print("a =", a);
	nvector_set_zero(a);
	nvector_print("a =", a);
	
	printf("\nmain: testing nvector_equal ...\n");
	nvector *c =nvector_alloc(2);

	nvector_set(c, 0, (double) 5.0);
	nvector_set(c, 1, (double) 1.0);

	nvector_print("a =", a);
	nvector_print("b =", b);
	nvector_print("c =", c);

	int a_b_equal = nvector_equal(a,b);
	int b_c_equal = nvector_equal(b,c);

	printf("a==b is %i\n", a_b_equal);
	printf("b==c is %i\n", b_c_equal);

	printf("\nmain: testing nvector_add ...\n");
	nvector *d = nvector_alloc(n);
	nvector *e = nvector_alloc(n);
	nvector *f = nvector_alloc(n);
	for (int i = 0; i < n; i++) {
		double x = RND, y = RND;
		nvector_set(d, i, x);
		nvector_set(e, i, y);
		nvector_set(f, i, x + y);
	}
	nvector_print("d =", d);
	nvector_print("e =", e);
	nvector_print("f =", f);
	nvector_add(d, e);
	nvector_print("d+e should   = ", f);
	nvector_print("d+e actually = ", d);

	if (nvector_equal(f, d)) printf("test passed\n");
	else printf("test failed\n");

	printf("\nmain: testing nvector_sub ...\n");
	nvector *d1 = nvector_alloc(n);
	nvector *e1 = nvector_alloc(n);
	nvector *f1 = nvector_alloc(n);
	for (int i = 0; i < n; i++) {
		double x1 = RND, y1 = RND;
		nvector_set(d1, i, x1);
		nvector_set(e1, i, y1);
		nvector_set(f1, i, x1 - y1);
	}
	nvector_print("d1 =", d1);
	nvector_print("e1 =", e1);
	nvector_print("f1 =", f1);
	nvector_sub(d1, e1);
	nvector_print("d1-e1 should   = ", f1);
	nvector_print("d1-e1 actually = ", d1);

	if (nvector_equal(f1, d1)) printf("test passed\n");
	else printf("test failed\n");

	printf("\nmain: testing nvector_scale ...\n");
	nvector_print("d1 =", d1);
	nvector_scale(d1, (double) 2.0);
	printf("x2 \n");
	nvector_print("d1 =", d1);

	nvector_free(v);
	nvector_free(a);
	nvector_free(b);
	nvector_free(c);
	nvector_free(d);
	nvector_free(e);
	nvector_free(f);
	nvector_free(d1);
	nvector_free(e1);
	nvector_free(f1);

	return 0;
}