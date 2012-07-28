#include "../real.h"
real_t chepolsum(int n, real_t x, real_t a[])
{
	int k;
	real_t h,r,s,tx;

	if (n == 0) return (a[0]);
	if (n == 1) return (a[0]+a[1]*x);
	tx=x+x;
	r=a[n];
	h=a[n-1]+r*tx;
	for (k=n-2; k>=1; k--) {
		s=r;
		r=h;
		h=a[k]+r*tx-s;
	}
	return (a[0]-r+h*x);
}
