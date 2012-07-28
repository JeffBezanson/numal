#include "../real.h"


void decsoltripiv(real_t sub[], real_t diag[], real_t super[], int n,
						real_t aux[], real_t b[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	int i,i1,n1,n2,*piv;
	real_t d,r,s,u,t,q,v,w,norm,norm1,norm2,tol,bi,bi1,bi2;

	piv=allocate_integer_vector(1,n);
	tol=aux[2];
	d=diag[1];
	r=super[1];
	bi=b[1];
	norm=norm2=fabs(d)+fabs(r);
	n2=n-2;
	for (i=1; i<=n2; i++) {
		i1=i+1;
		s=sub[i];
		t=diag[i1];
		q=super[i1];
		bi1=b[i1];
		norm1=norm2;
		norm2=fabs(s)+fabs(t)+fabs(q);
		if (norm2 > norm) norm=norm2;
		if (fabs(d)*norm2 < fabs(s)*norm1) {
			if (fabs(s) <= tol*norm2) {
				aux[3]=i-1;
				aux[5]=s;
				free_integer_vector(piv,1);
				return;
			}
			u=super[i]=t/s;
			b[i] = bi1 /= s;
			bi -= bi1*d;
			v=sub[i]=q/s;
			w = super[i1] = -v*d;
			d=diag[i1]=r-u*d;
			r=w;
			norm2=norm1;
			piv[i]=1;
		} else {
			if (fabs(d) <= tol*norm1) {
				aux[3]=i-1;
				aux[5]=d;
				free_integer_vector(piv,1);
				return;
			}
			u=super[i]=r/d;
			b[i] = bi /= d;
			bi=bi1-bi*s;
			d=diag[i1]=t-u*s;
			piv[i]=0;
			r=q;
		}
	}
	n1=n-1;
	s=sub[n1];
	t=diag[n];
	norm1=norm2;
	bi1=b[n];
	norm2=fabs(s)+fabs(t);
	if (norm2 > norm) norm=norm2;
	if (fabs(d)*norm2 < fabs(s)*norm1) {
		if (fabs(s) <= tol*norm2) {
			aux[3]=n2;
			aux[5]=s;
			free_integer_vector(piv,1);
			return;
		}
		u=super[n1]=t/s;
		b[n1] = bi1 /= s;
		bi -= bi1*d;
		d=r-u*d;
		norm2=norm1;
	} else {
		if (fabs(d) <= tol*norm1) {
			aux[3]=n2;
			aux[5]=d;
			free_integer_vector(piv,1);
			return;
		}
		u=super[n1]=r/d;
		b[n1] = bi /= d;
		bi=bi1-bi*s;
		d=t-u*s;
	}
	if (fabs(d) <= tol*norm2) {
		aux[3]=n1;
		aux[5]=d;
		free_integer_vector(piv,1);
		return;
	}
	aux[3]=n;
	aux[5]=norm;
	bi1=b[n]=bi/d;
	bi = b[n1] -= super[n1]*bi1;
	for (i=n-2; i>=1; i--) {
		bi2=bi1;
		bi1=bi;
		bi = b[i] -= super[i]*bi1 + ((piv[i]) ? sub[i]*bi2 : 0.0);
	}
	free_integer_vector(piv,1);
}
