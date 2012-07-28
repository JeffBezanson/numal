#include "../real.h"
void gsswts(int n, real_t zer[], real_t b[], real_t c[], real_t w[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void allortpol(int, real_t, real_t [], real_t [], real_t []);
	int j,k;
	real_t s,*p;

	p=allocate_real_vector(0,n-1);
	for (j=1; j<=n; j++) {
		allortpol(n-1,zer[j],b,c,p);
		s=0.0;
		for (k=n-1; k>=1; k--) s=(s+p[k]*p[k])/c[k];
		w[j]=1.0/(1.0+s);
	}
	free_real_vector(p,0);
}
