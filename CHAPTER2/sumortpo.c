#include "../real.h"
real_t sumortpol(int n, real_t x, real_t b[], real_t c[], real_t a[])
{
	int k;
	real_t h,r,s;

	if (n == 0) return (a[0]);
	r=a[n];
	s=0.0;
	for (k=n-1; k>=1; k--) {
		h=r;
		r=a[k]+(x-b[k])*r+s;
		s = -c[k]*h;
	}
	return (a[0]+(x-b[0])*r+s);
}
