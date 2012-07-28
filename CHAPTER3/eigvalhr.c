#include "../real.h"
void eigvalhrm(real_t **a, int n, int numval, real_t val[], real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshhrmtrival(real_t **, int, real_t [], real_t [], real_t []);
	void valsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t []);
	real_t *d,*bb;

	d=allocate_real_vector(1,n);
	bb=allocate_real_vector(1,n-1);
	hshhrmtrival(a,n,d,bb,em);
	valsymtri(d,bb,n,1,numval,val,em);
	free_real_vector(d,1);
	free_real_vector(bb,1);
}

