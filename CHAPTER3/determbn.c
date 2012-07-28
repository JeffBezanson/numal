#include "../real.h"


real_t determbnd(real_t a[], int n, int lw, int rw, int sgndet)
{
	int i,l;
	real_t p;

	l=1;
	p=1.0;
	lw += rw+1;
	for (i=1; i<=n; i++) {
		p=a[l]*p;
		l += lw;
	}
	return (fabs(p)*sgndet);
}
