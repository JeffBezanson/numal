#include "../real.h"


int maxelmrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x)
{
	int k;
	real_t r, s;

	s=0.0;
	for (k=l; k<=u; k++) {
		r=(a[i][k] += b[j][k]*x);
		if (fabs(r) > s) {
			s=fabs(r);
			l=k;
		}
	}
	return (l);
}
