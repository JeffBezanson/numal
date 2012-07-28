#include "../real.h"


void orthog(int n, int lc, int uc, real_t **x)
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	int i,j,k;
	real_t normx;

	for (j=lc; j<=uc; j++) {
		normx=sqrt(tammat(1,n,j,j,x,x));
		for (i=1; i<=n; i++) x[i][j] /=normx;
		for (k=j+1; k<=uc; k++) elmcol(1,n,k,j,x,x,-tammat(1,n,k,j,x,x));
	}
}

