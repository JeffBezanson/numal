#include "../real.h"
int qrivalhrm(real_t **a, int n, real_t val[], real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshhrmtrival(real_t **, int, real_t [], real_t [], real_t []);
	int qrivalsymtri(real_t [], real_t [], int, real_t []);
	int i;
	real_t *bb;

	bb=allocate_real_vector(1,n);
	hshhrmtrival(a,n,val,bb,em);
	bb[n]=0.0;
	i=qrivalsymtri(val,bb,n,em);
	free_real_vector(bb,1);
	return i;
}

