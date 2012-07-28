#include "../real.h"
real_t matvec(int l, int u, int i, real_t **a, real_t b[])
{
	int k;
	real_t s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[i][k]*b[k];
	return (s);
}
