#include "../real.h"
void chspol(int n, real_t a[])
{
	int k,l,twopow;

	if (n > 1) {
		for (k=0; k<=n-2; k++) {
			for (l=n-2; l>=k; l--) a[l] -= a[l+2];
			a[k+1] /= 2.0;
		}
		twopow=2;
		for (k=1; k<=n-2; k++) {
			a[k] *= twopow;
			twopow *= 2;
		}
		a[n-1] *= twopow;
		a[n] *= twopow;
	}
}
