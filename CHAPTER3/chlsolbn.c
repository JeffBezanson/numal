#include "../real.h"
void chlsolbnd(real_t a[], int n, int w, real_t b[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t scaprd1(int, int, int, int, int, real_t [], real_t []);
	int k,imax,kk,w1;

	kk = -w;
	w1=w+1;
	for (k=1; k<=n; k++) {
		kk += w1;
		b[k]=(b[k]-vecvec(((k <= w1) ? 1 : k-w),k-1,kk-k,b,a))/a[kk];
	}
	imax = -1;
	for (k=n; k>=1; k--) {
		if (imax < w) imax++;
		b[k]=(b[k]-scaprd1(kk+w,w,k+1,1,imax,a,b))/a[kk];
		kk -= w1;
	}
}
