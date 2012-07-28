#include "../real.h"
void sol(real_t **a, int n, int p[], real_t b[])
{
	real_t matvec(int, int, int, real_t **, real_t []);
	int k,pk;
	real_t r;

	for (k=1; k<=n; k++) {
		r=b[k];
		pk=p[k];
		b[k]=(b[pk]-matvec(1,k-1,k,a,b))/a[k][k];
		if (pk != k) b[pk]=r;
	}
	for (k=n; k>=1; k--) b[k] -= matvec(k+1,n,k,a,b);
}
