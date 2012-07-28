#include "../real.h"


void tfmsymtri2(real_t **a, int n, real_t d[], real_t b[], real_t bb[],
					real_t em[])
{
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void elmveccol(int, int, int, real_t [], real_t **, real_t);
	real_t tamvec(int, int, int, real_t **, real_t []);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	void elmcolvec(int, int, int, real_t **, real_t [], real_t);
	int i,j,r,r1;
	real_t w,x,a1,b0,bb0,machtol,norm;

	norm=0.0;
	for (j=1; j<=n; j++) {
		w=0.0;
		for (i=1; i<=j; i++) w += fabs(a[i][j]);
		for (i=j+1; i<=n; i++) w += fabs(a[j][i]);
		if (w > norm) norm=w;
	}
	machtol=em[0]*norm;
	em[1]=norm;
	r=n;
	for (r1=n-1; r1>=1; r1--) {
		d[r]=a[r][r];
		x=tammat(1,r-2,r,r,a,a);
		a1=a[r1][r];
		if (sqrt(x) <= machtol) {
			b0=b[r1]=a1;
			bb[r1]=b0*b0;
			a[r][r]=1.0;
		} else {
			bb0=bb[r1]=a1*a1+x;
			b0 = (a1 > 0.0) ? -sqrt(bb0) : sqrt(bb0);
			a1=a[r1][r]=a1-b0;
			w=a[r][r]=1.0/(a1*b0);
			for (j=1; j<=r1; j++)
				b[j]=(tammat(1,j,j,r,a,a)+matmat(j+1,r1,j,r,a,a))*w;
			elmveccol(1,r1,r,b,a,tamvec(1,r1,r,a,b)*w*0.5);
			for (j=1; j<=r1; j++) {
				elmcol(1,j,j,r,a,a,b[j]);
				elmcolvec(1,j,j,a,b,a[j][r]);
			}
			b[r1]=b0;
		}
		r=r1;
	}
	d[1]=a[1][1];
	a[1][1]=1.0;
	b[n]=bb[n]=0.0;
}
