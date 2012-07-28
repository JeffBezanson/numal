#include "../real.h"
void solbnd(real_t a[], int n, int lw, int rw, real_t m[],
				int p[], real_t b[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int f,i,k,kk,w,w1,w2,shift;
	real_t s;

	f=lw;
	shift = -lw;
	w1=lw-1;
	for (k=1; k<=n; k++) {
		if (f < n) f++;
		shift += w1;
		i=p[k];
		s=b[i];
		if (i != k) {
			b[i]=b[k];
			b[k]=s;
		}
		elmvec(k+1,f,shift,b,m,-s);
	}
	w1=lw+rw;
	w=w1+1;
	kk=(n+1)*w-w1;
	w2 = -1;
	shift=n*w1;
	for (k=n; k>=1; k--) {
		kk -= w;
		shift -= w1;
		if (w2 < w1) w2++;
		b[k]=(b[k]-vecvec(k+1,k+w2,shift,b,a))/a[kk];
	}
}
