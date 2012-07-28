#include "../real.h"
void lintfmpol(real_t p, real_t q, int n, real_t a[])
{
	void norderpol(int, int, real_t, real_t []);
	int k;
	real_t ppower;

	norderpol(n,n,q,a);
	ppower=p;
	for (k=1; k<=n; k++) {
		a[k] *= ppower;
		ppower *= p;
	}
}
