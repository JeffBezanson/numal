#include "../real.h"
void decinv(real_t **a, int n, real_t aux[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void dec(real_t **, int, real_t [], int []);
	void inv(real_t **, int, int []);
	int *p;

	p=allocate_integer_vector(1,n);
	dec(a,n,aux,p);
	if (aux[3] == n) inv(a,n,p);
	free_integer_vector(p,1);
}
