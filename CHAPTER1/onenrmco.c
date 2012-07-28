#include "../real.h"


real_t onenrmcol(int l, int u, int j, real_t **a)
{
	real_t sum;

	sum=0.0;
	for (; l<=u; l++) sum += fabs(a[l][j]);
	return (sum);
}
