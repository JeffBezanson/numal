#include "../real.h"
real_t pol(int n, real_t x, real_t a[])
{
	real_t r;

	r=0.0;
	for (; n>=0; n--) r=r*x+a[n];
	return (r);
}
