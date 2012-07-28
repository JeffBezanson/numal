#include "../real.h"
void baksymtri2(real_t **a, int n, int n1, int n2, real_t **vec)
{
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	int j,k;
	real_t w;

	for (j=2; j<=n; j++) {
		w=a[j][j];
		if (w < 0.0)
			for (k=n1; k<=n2; k++)
				elmcol(1,j-1,k,j,vec,a,tammat(1,j-1,j,k,a,vec)*w);
	}
}
