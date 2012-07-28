#include "../real.h"


real_t infnrmcol(int l, int u, int j, int *k, real_t **a)
{
	real_t r, max;

	max=0.0;
	*k=l;
	for (; l<=u; l++) {
		r=fabs(a[l][j]);
		if (r > max) {
			max=r;
			*k=l;
		}
	}
	return (max);
}
