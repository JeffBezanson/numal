#include "../real.h"
void lsqdglinv(real_t **a, int m, real_t aid[], int ci[], real_t diag[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int j,k,cik;
	real_t w;

	for (k=1; k<=m; k++) {
		diag[k]=1.0/aid[k];
		for (j=k+1; j<=m; j++) diag[j] = -tamvec(k,j-1,j,a,diag)/aid[j];
		diag[k]=vecvec(k,m,0,diag,diag);
	}
	for (k=m; k>=1; k--) {
		cik=ci[k];
		if (cik != k) {
			w=diag[k];
			diag[k]=diag[cik];
			diag[cik]=w;
		}
	}
}
