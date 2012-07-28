#include "../real.h"


void hshhrmtrival(real_t **a, int n, real_t d[], real_t bb[], real_t em[])
{
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t mattam(int, int, int, int, real_t **, real_t **);
	void elmveccol(int, int, int, real_t [], real_t **, real_t);
	void elmcolvec(int, int, int, real_t **, real_t [], real_t);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	void elmrow(int, int, int, int, real_t **, real_t **, real_t);
	void elmvecrow(int, int, int, real_t [], real_t **, real_t);
	void elmrowvec(int, int, int, real_t **, real_t [], real_t);
	void elmrowcol(int, int, int, int, real_t **, real_t **, real_t);
	void elmcolrow(int, int, int, int, real_t **, real_t **, real_t);
	int i,j,j1,jm1,r,rm1;
	real_t nrm,w,tol2,x,ar,ai,h,t,q,ajr,arj,dj,bbj,mod2;

	nrm=0.0;
	for (i=1; i<=n; i++) {
		w=fabs(a[i][i]);
		for (j=i-1; j>=1; j--) w += fabs(a[i][j])+fabs(a[j][i]);
		for (j=i+1; j<=n; j++) w += fabs(a[i][j])+fabs(a[j][i]);
		if (w > nrm) nrm=w;
	}
	t=em[0]*nrm;
	tol2=t*t;
	em[1]=nrm;
	r=n;
	for (rm1=n-1; rm1>=1; rm1--) {
		x=tammat(1,r-2,r,r,a,a)+mattam(1,r-2,r,r,a,a);
		ar=a[rm1][r];
		ai = -a[r][rm1];
		d[r]=a[r][r];
		if (x < tol2)
			bb[rm1]=ar*ar+ai*ai;
		else {
			mod2=ar*ar+ai*ai;
			if (mod2 == 0.0) {
				a[rm1][r]=sqrt(x);
				t=x;
			} else {
				x += mod2;
				h=sqrt(mod2*x);
				t=x+h;
				h=1.0+x/h;
				a[r][rm1] = -ai*h;
				a[rm1][r]=ar*h;
			}
			j=1;
			jm1=0;
			for (j1=2; j1<=r; j1++) {
				d[j]=(tammat(1,j,j,r,a,a)+matmat(j1,rm1,j,r,a,a)+
						mattam(1,jm1,j,r,a,a)-matmat(j1,rm1,r,j,a,a))/t;
				bb[j]=(matmat(1,jm1,j,r,a,a)-tammat(j1,rm1,j,r,a,a)-
						matmat(1,j,r,j,a,a)-mattam(j1,rm1,j,r,a,a))/t;
				jm1=j;
				j=j1;
			}
			q=(tamvec(1,rm1,r,a,d)-matvec(1,rm1,r,a,bb))/t/2.0;
			elmveccol(1,rm1,r,d,a,-q);
			elmvecrow(1,rm1,r,bb,a,q);
			j=1;
			for (j1=2; j1<=r; j1++) {
				ajr=a[j][r];
				arj=a[r][j];
				dj=d[j];
				bbj=bb[j];
				elmrowvec(j,rm1,j,a,d,-ajr);
				elmrowvec(j,rm1,j,a,bb,arj);
				elmrowcol(j,rm1,j,r,a,a,-dj);
				elmrow(j,rm1,j,r,a,a,bbj);
				elmcolvec(j1,rm1,j,a,d,-arj);
				elmcolvec(j1,rm1,j,a,bb,-ajr);
				elmcol(j1,rm1,j,r,a,a,bbj);
				elmcolrow(j1,rm1,j,r,a,a,dj);
				j=j1;
			}
			bb[rm1]=x;
		}
		r=rm1;
	}
	d[1]=a[1][1];
}
