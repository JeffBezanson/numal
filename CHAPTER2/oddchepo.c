#include "../real.h"
real_t oddchepolsum(int n, real_t x, real_t a[])
{
	int k;
	real_t h,r,s,y;

	if (n == 0) return (x*a[0]);
	if (n == 1) return (x*(a[0]+a[1]*(4.0*x*x-3.0)));
	y=4.0*x*x-2.0;
	r=a[n];
	h=a[n-1]+r*y;
	for (k=n-2; k>=0; k--) {
		s=r;
		r=h;
		h=a[k]+r*y-s;
	}
	return (x*(h-r));
}
