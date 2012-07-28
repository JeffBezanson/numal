#include "../real.h"
void hestgl2(int n, real_t **a, real_t **b)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hsh2col(int, int, int, int, real_t, real_t, real_t **, real_t **);
	void hsh2row2(int, int, int, int, real_t, real_t,
						real_t **, real_t **);
	int nm1,k,l,k1,l1;

	if (n > 2) {
		for (k=2; k<=n; k++)
			for (l=1; l<=k-1; l++) b[k][l]=0.0;
		nm1=n-1;
		k=1;
		for (k1=2; k1<=nm1; k1++) {
			l1=n;
			for (l=n-1; l>=k1; l--) {
				hsh2col(k,l,n,l,a[l][k],a[l1][k],a,b);
				a[l1][k]=0.0;
				hsh2row2(1,n,l1,l,b[l1][l1],b[l1][l],a,b);
				b[l1][l]=0.0;
				l1=l;
			}
			k=k1;
		}
	}
}
