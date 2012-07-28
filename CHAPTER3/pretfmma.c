#include "../real.h"
void pretfmmat(real_t **a, int m, int n, real_t d[])
{
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	int i,i1,j;
	real_t g,h;

	for (i=n; i>=1; i--) {
		i1=i+1;
		g=d[i];
		h=g*a[i][i];
		for (j=i1; j<=n; j++) a[i][j]=0.0;
		if (h < 0.0) {
			for (j=i1; j<=n; j++)
				elmcol(i,m,j,i,a,a,tammat(i1,m,i,j,a,a)/h);
			for (j=i; j<=m; j++) a[j][i] /= g;
		} else
			for (j=i; j<=m; j++) a[j][i]=0.0;
		a[i][i] += 1.0;
	}
}
