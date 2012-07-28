#include "../real.h"


void eialpha(real_t x, int n, real_t alpha[])
{
	int k;
	real_t a,b,c;

	c=1.0/x;
	a=exp(-x);
	b=alpha[0]=a*c;
	for (k=1; k<=n; k++) alpha[k]=b=(a+k*b)*c;
}
