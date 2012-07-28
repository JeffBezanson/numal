#include "../real.h"
real_t tamvec(int l, int u, int i, real_t **a, real_t b[])
{
	int k;
	real_t s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[k][i]*b[k];
	return (s);
}
