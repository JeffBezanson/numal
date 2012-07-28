#include "../real.h"


void dectri(real_t sub[], real_t diag[], real_t super[],
				int n, real_t aux[])
{
	int i,n1;
	real_t d,r,s,u,norm,norm1,tol;

	tol=aux[2];
	d=diag[1];
	r=super[1];
	norm=norm1=fabs(d)+fabs(r);
	if (fabs(d) <= norm1*tol) {
		aux[3]=0.0;
		aux[5]=d;
		return;
	}
	u=super[1]=r/d;
	s=sub[1];
	n1=n-1;
	for (i=2; i<=n1; i++) {
		d=diag[i];
		r=super[i];
		norm1=fabs(s)+fabs(d)+fabs(r);
		diag[i] = d -= u*s;
		if (fabs(d) <= norm1*tol) {
			aux[3]=i-1;
			aux[5]=d;
			return;
		}
		u=super[i]=r/d;
		s=sub[i];
		if (norm1 > norm) norm=norm1;
	}
	d=diag[n];
	norm1=fabs(d)+fabs(s);
	diag[n] = d -= u*s;
	if (fabs(d) <= norm1*tol) {
		aux[3]=n1;
		aux[5]=d;
		return;
	}
	if (norm1 > norm) norm=norm1;
	aux[3]=n;
	aux[5]=norm;
}
