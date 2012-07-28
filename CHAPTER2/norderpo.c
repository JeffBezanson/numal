#include "../real.h"
void norderpol(int n, int k, real_t x, real_t a[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int i,j,nm1;
	real_t xj,aa,h,*xx;

	if (x != 0.0) {
		xx=allocate_real_vector(0,n);
		xj=1;
		for (j=1; j<=n; j++) {
			xx[j] = xj *= x;
			a[j] *= xj;
		}
		h=aa=a[n];
		nm1=n-1;
		for (i=nm1; i>=0; i--) h = a[i] += h;
		for (j=1; j<=k; j++) {
			h=aa;
			for (i=nm1; i>=j; i--) h = a[i] += h;
			a[j]=h/xx[j];
		}
		free_real_vector(xx,0);
	}
}
