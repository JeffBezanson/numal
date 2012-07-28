#include "../real.h"
void chlsol2(real_t **a, int n, real_t b[])
{
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int i;

	for (i=1; i<=n; i++) b[i]=(b[i]-tamvec(1,i-1,i,a,b))/a[i][i];
	for (i=n; i>=1; i--) b[i]=(b[i]-matvec(i+1,n,i,a,b))/a[i][i];
}
