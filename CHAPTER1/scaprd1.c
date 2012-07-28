#include "../real.h"
real_t scaprd1(int la, int sa, int lb, int sb, int n, real_t a[], real_t b[])
{
	int k;
	real_t s;

	s=0.0;
	for (k=1; k<=n; k++) {
		s += a[la]*b[lb];
		la += sa;
		lb += sb;
	}
	return (s);
}
