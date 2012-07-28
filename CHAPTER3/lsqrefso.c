#include "../real.h"


void lsqrefsol(real_t **a, real_t **qr, int n, int m, int n1,
					real_t aux[], real_t aid[], int ci[], real_t b[],
					real_t *ldx, real_t x[], real_t res[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	void ichcol(int, int, int, int, real_t **);
	void elmveccol(int, int, int, real_t [], real_t **, real_t);
	int i,j,k,s,startup;
	real_t c1,nexve,ndx,ndr,d,corrnorm,*f,*g;
	double dtemp;

	f=allocate_real_vector(1,n);
	g=allocate_real_vector(1,m);
	for (j=1; j<=m; j++) {
		s=ci[j];
		if (s != j) ichcol(1,n,j,s,a);
	}
	for (j=1; j<=m; j++) x[j]=g[j]=0.0;
	for (i=1; i<=n; i++) {
		res[i]=0.0;
		f[i]=b[i];
	}
	k=0;
	do {
		startup = (k <= 1);
		ndx=ndr=0.0;
		if (k != 0) {
			for (i=1; i<=n; i++) res[i] += f[i];
			for (s=1; s<=m; s++) {
				x[s] += g[s];
				dtemp=0.0;
				for (i=1; i<=n; i++)
					dtemp += (double)a[i][s]*(double)res[i];
				d=dtemp;
				g[s]=(-d-tamvec(1,s-1,s,qr,g))/aid[s];
			}
			for (i=1; i<=n; i++) {
				dtemp = (i > n1) ? res[i] : 0.0;
				for (s=1; s<=m; s++)
					dtemp += (double)a[i][s]*(double)x[s];
				f[i]=(double)b[i]-dtemp;
			}
		}
		nexve=sqrt(vecvec(1,m,0,x,x)+vecvec(1,n,0,res,res));
		for (s=1; s<=n1; s++)
			elmveccol(s,n1,s,f,qr,tamvec(s,n1,s,qr,f)/(qr[s][s]*aid[s]));
		for (i=n1+1; i<=n; i++)
			f[i] -= matvec(1,n1,i,qr,f);
		for (s=n1+1; s<=m; s++)
			elmveccol(s,n,s,f,qr,tamvec(s,n,s,qr,f)/(qr[s][s]*aid[s]));
		for (i=1; i<=m; i++) {
			c1=f[i];
			f[i]=g[i];
			g[i] = (i > n1) ? c1-g[i] : c1;
		}
		for (s=m; s>=1; s--) {
			g[s]=(g[s]-matvec(s+1,m,s,qr,g))/aid[s];
			ndx += g[s]*g[s];
		}
		for (s=m; s>=n1+1; s--)
			elmveccol(s,n,s,f,qr,tamvec(s,n,s,qr,f)/(qr[s][s]*aid[s]));
		for (s=1; s<=n1; s++)
			f[s] -= tamvec(n1+1,n,s,qr,f);
		for (s=n1; s>=1; s--)
			elmveccol(s,n1,s,f,qr,tamvec(s,n1,s,qr,f)/(qr[s][s]*aid[s]));
		aux[7]=k;
		for (i=1; i<=n; i++) ndr += f[i]*f[i];
		corrnorm=sqrt(ndx+ndr);
		k++;
	} while (startup || (corrnorm>aux[2]*nexve && k<=aux[6]));
	*ldx=sqrt(ndx);
	for (s=m; s>=1; s--) {
		j=ci[s];
		if (j != s) {
			c1=x[j];
			x[j]=x[s];
			x[s]=c1;
			ichcol(1,n,j,s,a);
		}
	}
	free_real_vector(f,1);
	free_real_vector(g,1);
}
