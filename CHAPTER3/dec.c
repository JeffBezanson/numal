#include "../real.h"


void dec(real_t **a, int n, real_t aux[], int p[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t mattam(int, int, int, int, real_t **, real_t **);
	void ichrow(int, int, int, int, real_t **);
	int i,k,k1,pk,d;
	real_t r,s,eps,*v;

	v=allocate_real_vector(1,n);
	r = -1.0;
	for (i=1; i<=n; i++) {
		s=sqrt(mattam(1,n,i,i,a,a));
		if (s > r) r=s;
		v[i]=1.0/s;
	}
	eps=aux[2]*r;
	d=1.0;
	for (k=1; k<=n; k++) {
		r = -1.0;
		k1=k-1;
		for (i=k; i<=n; i++) {
			a[i][k] -= matmat(1,k1,i,k,a,a);
			s=fabs(a[i][k])*v[i];
			if (s > r) {
				r=s;
				pk=i;
			}
		}
		p[k]=pk;
		v[pk]=v[k];
		s=a[pk][k];
		if (fabs(s) < eps) break;
		if (s < 0.0) d = -d;
		if (pk != k) {
			d = -d;
			ichrow(1,n,k,pk,a);
		}
		for (i=k+1; i<=n; i++) a[k][i]=(a[k][i]-matmat(1,k1,k,i,a,a))/s;
	}
	aux[1]=d;
	aux[3]=k-1;
	free_real_vector(v,1);
}
