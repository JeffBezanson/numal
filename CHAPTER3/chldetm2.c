#include "../real.h"
real_t chldeterm2(real_t **a, int n)
{
	int k;
	real_t d;

	d=1.0;
	for (k=1; k<=n; k++) d *= a[k][k];
	return (d*d);
}
