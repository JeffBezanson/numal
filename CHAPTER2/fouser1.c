#include "../real.h"


real_t fouser1(int n, real_t theta, real_t a[], real_t b[])
{
	int i;
	real_t r,s,h,co,si;

	r=s=0.0;
	co=cos(theta);
	si=sin(theta);
	for (i=n; i>=1; i--) {
		h=co*r+si*s+a[i];
		s=co*s-si*r+b[i];
		r=h;
	}
	return (co*r+si*s+a[0]);
}
