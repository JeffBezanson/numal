#include "../real.h"
real_t ortpol(int n, real_t x, real_t b[], real_t c[])
{
	int k,l;
	real_t r,s,h;

	if (n == 0) return (1.0);
	r=x-b[0];
	s=1.0;
	l=n-1;
	for (k=1; k<=l; k++) {
		h=r;
		r=(x-b[k])*r-c[k]*s;
		s=h;
	}
	return (r);
}
