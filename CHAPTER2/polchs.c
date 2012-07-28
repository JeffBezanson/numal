#include "../real.h"
void polchs(int n, real_t a[])
{
	int k,l,twopow;

	if (n > 1) {
		twopow=2;
		for (k=1; k<=n-2; k++) {
			a[k] /= twopow;
			twopow *= 2;
		}
		a[n-1]=2.0*a[n-1]/twopow;
		a[n] /= twopow;
		a[n-2] += a[n];
		for (k=n-2; k>=1; k--) {
			a[k-1] += a[k+1];
			a[k]=2.0*a[k]+a[k+2];
			for (l=k+1; l<=n-2; l++) a[l] += a[l+2];
		}
	}
}
