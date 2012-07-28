#include "../real.h"


void dectripiv(real_t sub[], real_t diag[], real_t super[], int n,
					real_t aid[], real_t aux[], int piv[])
{
	int i,i1,n1,n2;
	real_t d,r,s,u,t,q,v,w,norm,norm1,norm2,tol;

	tol=aux[2];
	d=diag[1];
	r=super[1];
	norm=norm2=fabs(d)+fabs(r);
	n2=n-2;
	for (i=1; i<=n2; i++) {
		i1=i+1;
		s=sub[i];
		t=diag[i1];
		q=super[i1];
		norm1=norm2;
		norm2=fabs(s)+fabs(t)+fabs(q);
		if (norm2 > norm) norm=norm2;
		if (fabs(d)*norm2 < fabs(s)*norm1) {
			if (fabs(s) <= tol*norm2) {
				aux[3]=i-1;
				aux[5]=s;
				return;
			}
			diag[i]=s;
			u=super[i]=t/s;
			v=aid[i]=q/s;
			sub[i]=d;
			w = super[i1] = -v*d;
			d=diag[i1]=r-u*d;
			r=w;
			norm2=norm1;
			piv[i]=1;
		} else {
			if (fabs(d) <= tol*norm1) {
				aux[3]=i-1;
				aux[5]=d;
				return;
			}
			u=super[i]=r/d;
			d=diag[i1]=t-u*s;
			aid[i]=0.0;
			piv[i]=0;
			r=q;
		}
	}
	n1=n-1;
	s=sub[n1];
	t=diag[n];
	norm1=norm2;
	norm2=fabs(s)+fabs(t);
	if (norm2 > norm) norm=norm2;
	if (fabs(d)*norm2 < fabs(s)*norm1) {
		if (fabs(s) <= tol*norm2) {
			aux[3]=n2;
			aux[5]=s;
			return;
		}
		diag[n1]=s;
		u=super[n1]=t/s;
		sub[n1]=d;
		d=diag[n]=r-u*d;
		norm2=norm1;
		piv[n1]=1;
	} else {
		if (fabs(d) <= tol*norm1) {
			aux[3]=n2;
			aux[5]=d;
			return;
		}
		u=super[n1]=r/d;
		d=diag[n]=t-u*s;
		piv[n1]=0;
	}
	if (fabs(d) <= tol*norm2) {
		aux[3]=n1;
		aux[5]=d;
		return;
	}
	aux[3]=n;
	aux[5]=norm;
}
