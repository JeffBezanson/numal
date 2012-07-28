#include "../real.h"


void gsselm(real_t **a, int n, real_t aux[], int ri[], int ci[])
{
	void elmrow(int, int, int, int, real_t **, real_t **, real_t);
	int maxelmrow(int, int, int, int, real_t **, real_t **, real_t);
	void ichcol(int, int, int, int, real_t **);
	void ichrow(int, int, int, int, real_t **);
	void rowcst(int, int, int, real_t **, real_t);
	real_t absmaxmat(int, int, int, int, int *, int *, real_t **);
	int i,j,p,q,r,r1,jpiv,rank,signdet,partial;
	real_t crit,pivot,rgrow,max,aid,max1,eps;

	aux[5]=rgrow=absmaxmat(1,n,1,n,&i,&j,a);
	crit=n*rgrow*aux[4];
	eps=rgrow*aux[2];
	max=0.0;
	rank=n;
	signdet=1;
	partial = rgrow != 0;
	for (q=1; q<=n; q++)
		if (q != j) {
			aid=fabs(a[i][q]);
			if (aid > max) max=aid;
		}
	rgrow += max;
	for (r=1; r<=n; r++) {
		r1=r+1;
		if (i != r) {
			signdet = -signdet;
			ichrow(1,n,r,i,a);
		}
		if (j != r) {
			signdet = -signdet;
			ichcol(1,n,r,j,a);
		}
		ri[r]=i;
		ci[r]=j;
		pivot=a[r][r];
		if (pivot < 0.0) signdet = -signdet;
		if (partial) {
			max=max1=0.0;
			j=r1;
			rowcst(r1,n,r,a,1.0/pivot);
			for (p=r1; p<=n; p++) {
				elmrow(r1,n,p,r,a,a,-a[p][r]);
				aid=fabs(a[p][r1]);
				if (max < aid) {
					max=aid;
					i=p;
				}
			}
			for (q=r1+1; q<=n; q++) {
				aid=fabs(a[i][q]);
				if (max1 < aid) max1=aid;
			}
			aid=rgrow;
			rgrow += max1;
			if ((rgrow > crit) || (max < eps)) {
				partial=0;
				rgrow=aid;
				max=absmaxmat(r1,n,r1,n,&i,&j,a);
			}
		} else {
			if (max <= eps) {
				rank=r-1;
				if (pivot < 0.0) signdet = -signdet;
				break;
			}
			max = -1.0;
			rowcst(r1,n,r,a,1.0/pivot);
			for (p=r1; p<=n; p++) {
				jpiv=maxelmrow(r1,n,p,r,a,a,-a[p][r]);
				aid=fabs(a[p][jpiv]);
				if (max < aid) {
					max=aid;
					i=p;
					j=jpiv;
				}
			}
			if (rgrow < max) rgrow=max;
		}
	}
	aux[1]=signdet;
	aux[3]=rank;
	aux[7]=rgrow;
}
