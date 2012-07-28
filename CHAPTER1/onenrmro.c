#include "../real.h"


real_t onenrmrow(int l, int u, int i, real_t **a)
{
	real_t sum;

	sum=0.0;
	for (; l<=u; l++) sum += fabs(a[i][l]);
	return (sum);
}
