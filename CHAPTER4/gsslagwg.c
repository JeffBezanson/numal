#include "../real.h"


void gsslagwghts(int n, real_t alfa, real_t x[], real_t w[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void alllagzer(int, real_t, real_t []);
	real_t gamma(real_t);
	int i,j;
	real_t h0,s,r0,r1,r2,xi,*a,*b;

	a=allocate_real_vector(0,n);
	b=allocate_real_vector(0,n);
	a[0]=1.0+alfa;
	a[1]=3.0+alfa;
	b[1]=sqrt(a[0]);
	for (i=2; i<=n-1; i++) {
		a[i]=i+i+alfa+1.0;
		b[i]=sqrt(i*(i+alfa));
	}
	alllagzer(n,alfa,x);
	h0=gamma(1.0+alfa);
	for (i=1; i<=n; i++) {
		xi=x[i];
		r0=1.0;
		r1=(xi-a[0])/b[1];
		s=1.0+r1*r1;
		for (j=2; j<=n-1; j++) {
			r2=((xi-a[j-1])*r1-b[j-1]*r0)/b[j];
			r0=r1;
			r1=r2;
			s += r2*r2;
		}
		w[i]=h0/s;
	}
	free_real_vector(a,0);
	free_real_vector(b,0);
}
