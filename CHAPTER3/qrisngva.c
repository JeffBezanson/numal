#include "../real.h"
int qrisngval(real_t **a, int m, int n, real_t val[], real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshreabid(real_t **, int, int, real_t [], real_t [], real_t []);
	int qrisngvalbid(real_t [], real_t [], int, real_t []);
	int i;
	real_t *b;

	b=allocate_real_vector(1,n);
	hshreabid(a,m,n,val,b,em);
	i=qrisngvalbid(val,b,n,em);
	free_real_vector(b,1);
	return i;
}
