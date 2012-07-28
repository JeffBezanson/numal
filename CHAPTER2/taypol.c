#include "../real.h"
void taypol(int n, int k, real_t x, real_t a[])
{
	int i,j,nm1;
	real_t xj,aa,h;

	if (x != 0.0) {
		xj=1;
		for (j=1; j<=n; j++) {
			xj *= x;
			a[j] *= xj;
		}
		aa=a[n];
		nm1=n-1;
		for (j=0; j<=k; j++) {
			h=aa;
			for (i=nm1; i>=j; i--) h = a[i] += h;
		}
	} else {
		for (; k>=1; n--) a[k]=0;
	}
}
