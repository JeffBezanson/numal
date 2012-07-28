#include "../real.h"
real_t vecvec(int l, int u, int shift, real_t a[], real_t b[])
{
	int k;
	real_t s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[k]*b[k+shift];
	return (s);
}
