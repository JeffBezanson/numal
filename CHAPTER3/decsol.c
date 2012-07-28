#include "../real.h"
void decsol(real_t **a, int n, real_t aux[], real_t b[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void sol(real_t **, int, int [], real_t []);
	void dec(real_t **, int, real_t [], int []);
	int *p;

	p=allocate_integer_vector(1,n);
	dec(a,n,aux,p);
	if (aux[3] == n) sol(a,n,p,b);
	free_integer_vector(p,1);
}
