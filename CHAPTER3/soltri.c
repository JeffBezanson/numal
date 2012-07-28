#include "../real.h"
void soltri(real_t sub[], real_t diag[], real_t super[], int n, real_t b[])
{
	int i;
	real_t r;

	r = b[1] /= diag[1];
	for (i=2; i<=n; i++) r=b[i]=(b[i]-sub[i-1]*r)/diag[i];
	for (i=n-1; i>=1; i--) r = b[i] -= super[i]*r;
}
