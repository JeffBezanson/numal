#include "../real.h"
real_t jfrac(int n, real_t a[], real_t b[])
{
	int i;
	real_t d;

	d=0.0;
	for (i=n; i>=1; i--) d=a[i]/(b[i]+d);
	return (d+b[0]);
}
