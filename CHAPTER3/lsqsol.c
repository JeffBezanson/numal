#include "../real.h"
void lsqsol(real_t **a, int n, int m, real_t aid[], int ci[], real_t b[])
{
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	void elmveccol(int, int, int, real_t [], real_t **, real_t);
	int k,cik;
	real_t w;

	for (k=1; k<=m; k++)
		elmveccol(k,n,k,b,a,tamvec(k,n,k,a,b)/(aid[k]*a[k][k]));
	for (k=m; k>=1; k--) b[k]=(b[k]-matvec(k+1,m,k,a,b))/aid[k];
	for (k=m; k>=1; k--) {
		cik=ci[k];
		if (cik != k) {
			w=b[k];
			b[k]=b[cik];
			b[cik]=w;
		}
	}
}
