#include "../real.h"


real_t inv1(real_t **a, int n, int ri[], int ci[], int withnorm)
{
	void ichrow(int, int, int, int, real_t **);
	void inv(real_t **, int, int []);
	int l,k,k1;
	real_t aid,nrminv;

	inv(a,n,ri);
	nrminv=0.0;
	if (withnorm)
		for (l=1; l<=n; l++) nrminv += fabs(a[l][n]);
	for (k=n-1; k>=1; k--) {
		if (withnorm) {
			aid=0.0;
			for (l=1; l<=n; l++) aid += fabs(a[l][k]);
			if (nrminv < aid) nrminv=aid;
		}
		k1=ci[k];
		if (k1 != k) ichrow(1,n,k,k1,a);
	}
	return (nrminv);
}
