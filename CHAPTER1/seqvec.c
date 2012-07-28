#include "../real.h"
real_t seqvec(int l, int u, int il, int shift, real_t a[], real_t b[])
{
	real_t s;

	s=0.0;
	for (; l<=u; l++) {
		s += a[il]*b[l+shift];
		il += l;
	}
	return (s);
}
