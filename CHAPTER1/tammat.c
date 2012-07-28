#include "../real.h"
real_t tammat(int l, int u, int i, int j, real_t **a, real_t **b)
{
	int k;
	real_t s;

	s=0.0;
	for (k=l; k<=u; k++) s += a[k][i]*b[k][j];
	return (s);
}
