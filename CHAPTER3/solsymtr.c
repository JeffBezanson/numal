#include "../real.h"
void solsymtri(real_t diag[], real_t co[], int n, real_t b[])
{
	int i;
	real_t r,s;

	r=b[1];
	b[1]=r/diag[1];
	for (i=2; i<=n; i++) {
		r=b[i]-co[i-1]*r;
		b[i]=r/diag[i];
	}
	s=b[n];
	for (i=n-1; i>=1; i--) s = b[i] -= co[i]*s;
}
