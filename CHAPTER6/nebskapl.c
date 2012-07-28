#include "../real.h"


void nonexpbesskaplusn(real_t a, real_t x, int nmax, real_t kan[])
{
	void nonexpbesska01(real_t, real_t, real_t *, real_t *);
	int n;
	real_t k1;

	nonexpbesska01(a,x,&kan[0],&k1);
	a -= 1.0;
	x=2.0/x;
	if (nmax > 0) kan[1]=k1;
	for (n=2; n<=nmax; n++) kan[n]=kan[n-2]+(a+n)*x*kan[n-1];
}
