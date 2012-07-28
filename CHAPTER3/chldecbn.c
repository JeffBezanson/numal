#include "../real.h"


void chldecbnd(real_t a[], int n, int w, real_t aux[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	int j,k,jmax,kk,kj,w1,start;
	real_t r,eps,max;

	max=0.0;
	kk = -w;
	w1=w+1;
	for (j=1; j<=n; j++) {
		kk += w1;
		if (a[kk] > max) max=a[kk];
	}
	jmax=w;
	w1=w+1;
	kk = -w;
	eps=aux[2]*max;
	for (k=1; k<=n; k++) {
		if (k+w > n) jmax--;
		kk += w1;
		start=kk-k+1;
		r=a[kk]-vecvec(((k <= w1) ? start : kk-w),kk-1,0,a,a);
		if (r <= eps) {
			aux[3]=k-1;
			return;
		}
		a[kk]=r=sqrt(r);
		kj=kk;
		for (j=1; j<=jmax; j++) {
			kj += w;
			a[kj]=(a[kj]-vecvec(((k+j <= w1) ? start : kk-w+j),
							kk-1,kj-kk,a,a))/r;
		}
	}
	aux[3]=n;
}
