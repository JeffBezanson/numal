#include "../real.h"


void bessyaplusn(real_t a, real_t x, int nmax, real_t yan[])
{
	void bessya01(real_t, real_t, real_t *, real_t *);
	int n;
	real_t y1;

	bessya01(a,x,&yan[0],&y1);
	a -= 1.0;
	x=2.0/x;
	if (nmax > 0) yan[1]=y1;
	for (n=2; n<=nmax; n++) yan[n] = -yan[n-2]+(a+n)*x*yan[n-1];
}
