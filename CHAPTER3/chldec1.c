#include "../real.h"


void chldec1(real_t a[], int n, real_t aux[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	int j,k,kk,kj,low,up;
	real_t r,epsnorm;

	r=0.0;
	kk=0;
	for (k=1; k<=n; k++) {
		kk += k;
		if (a[kk] > r) r=a[kk];
	}
	epsnorm=aux[2]*r;
	kk=0;
	for (k=1; k<=n; k++) {
		kk += k;
		low=kk-k+1;
		up=kk-1;
		r=a[kk]-vecvec(low,up,0,a,a);
		if (r <= epsnorm) {
			aux[3]=k-1;
			return;
		}
		a[kk]=r=sqrt(r);
		kj=kk+k;
		for (j=k+1; j<=n; j++) {
			a[kj]=(a[kj]-vecvec(low,up,kj-kk,a,a))/r;
			kj +=j;
		}
	}
	aux[3]=n;
}
