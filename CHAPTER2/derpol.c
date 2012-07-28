#include "../real.h"
void derpol(int n, int k, real_t x, real_t a[])
{
	void norderpol(int, int, real_t, real_t []);
	int j;
	real_t fac;

	fac=1.0;
	norderpol(n,k,x,a);
	for (j=2; j<=k; j++) {
		fac *= j;
		a[j] *=fac;
	}
}
