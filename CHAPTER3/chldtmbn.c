#include "../real.h"
real_t chldetermbnd(real_t a[], int n, int w)
{
	int j,kk,w1;
	real_t p;

	w1=w+1;
	kk = -w;
	p=1.0;
	for (j=1; j<=n; j++) {
		kk += w1;
		p *= a[kk];
	}
	return (p*p);
}
