#include "../real.h"
int qrivalsym2(real_t **a, int n, real_t val[], real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void tfmsymtri2(real_t **, int, real_t [], real_t [], real_t [],
						real_t []);
	int qrivalsymtri(real_t [], real_t [], int, real_t []);
	int i;
	real_t *b,*bb;

	b=allocate_real_vector(1,n);
	bb=allocate_real_vector(1,n);
	tfmsymtri2(a,n,val,b,bb,em);
	i=qrivalsymtri(val,bb,n,em);
	free_real_vector(b,1);
	free_real_vector(bb,1);
	return i;
}

