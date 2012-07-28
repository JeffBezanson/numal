#include "../real.h"


real_t infnrmrow(int l, int u, int i, int *k, real_t **a)
{
	real_t r, max;

	max=0.0;
	*k=l;
	for (; l<=u; l++) {
		r=fabs(a[i][l]);
		if (r > max) {
			max=r;
			*k=l;
		}
	}
	return (max);
}
