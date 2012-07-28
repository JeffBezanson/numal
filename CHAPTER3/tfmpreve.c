#include "../real.h"
void tfmprevec(real_t **a, int n)
{
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	int i,j,j1,k;
	real_t ab;

	j1=1;
	for (j=2; j<=n; j++) {
		for (i=1; i<=j1-1; i++) a[i][j1]=0.0;
		for (i=j; i<=n; i++) a[i][j1]=0.0;
		a[j1][j1]=1.0;
		ab=a[j][j];
		if (ab < 0)
			for (k=1; k<=j1; k++)
				elmcol(1,j1,k,j,a,a,tammat(1,j1,j,k,a,a)*ab);
		j1=j;
	}
	for (i=n-1; i>=1; i--) a[i][n]=0.0;
	a[n][n]=1.0;
}
