#include "../real.h"


void decsymtri(real_t diag[], real_t co[], int n, real_t aux[])
{
	int i,n1;
	real_t d,r,s,u,tol,norm,normr;

	tol=aux[2];
	d=diag[1];
	r=co[1];
	norm=normr=fabs(d)+fabs(r);
	if (fabs(d) <= normr*tol) {
		aux[3]=0.0;
		aux[5]=d;
		return;
	}
	u=co[1]=r/d;
	n1=n-1;
	for (i=2; i<=n1; i++) {
		s=r;
		r=co[i];
		d=diag[i];
		normr=fabs(s)+fabs(d)+fabs(r);
		diag[i] = d -= u*s;
		if (fabs(d) <= normr*tol) {
			aux[3]=i-1;
			aux[5]=d;
			return;
		}
		u=co[i]=r/d;
		if (normr > norm) norm=normr;
	}
	d=diag[n];
	normr=fabs(d)+fabs(r);
	diag[n] = d -= u*s;
	if (fabs(d) <= normr*tol) {
		aux[3]=n1;
		aux[5]=d;
		return;
	}
	if (normr > norm) norm=normr;
	aux[3]=n;
	aux[5]=norm;
}
