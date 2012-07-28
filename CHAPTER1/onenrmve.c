#include "../real.h"


real_t onenrmvec(int l, int u, real_t a[])
{
	real_t sum;

	sum=0.0;
	for (; l<=u; l++) sum += fabs(a[l]);
	return (sum);
}
