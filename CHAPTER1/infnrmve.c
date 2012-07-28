#include "../real.h"


real_t infnrmvec(int l, int u, int *k, real_t a[])
{
	real_t r, max;

	max=0.0;
	*k=l;
	for (; l<=u; l++) {
		r=fabs(a[l]);
		if (r > max) {
			max=r;
			*k=l;
		}
	}
	return (max);
}
