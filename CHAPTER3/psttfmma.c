#include "../real.h"
void psttfmmat(real_t **a, int n, real_t **v, real_t b[])
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	int i,i1,j;
	real_t h;

	i1=n;
	v[n][n]=1.0;
	for (i=n-1; i>=1; i--) {
		h=b[i]*a[i][i1];
		if (h < 0.0) {
			for (j=i1; j<=n; j++) v[j][i]=a[i][j]/h;
			for (j=i1; j<=n; j++)
				elmcol(i1,n,j,i,v,v,matmat(i1,n,i,j,a,v));
		}
		for (j=i1; j<=n; j++) v[i][j]=v[j][i]=0.0;
		v[i][i]=1.0;
		i1=i;
	}
}
