#include "../real.h"

void alllagzer(int n, real_t alfa, real_t zer[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void allzerortpol(int, real_t [], real_t [], real_t [], real_t []);
	int i;
	real_t *a,*b,em[6];

	a=allocate_real_vector(0,n);
	b=allocate_real_vector(0,n);
	b[0]=0.0;
	a[n-1]=n+n+alfa-1.0;
	for (i=1; i<=n-1; i++) {
		a[i-1]=i+i+alfa-1.0;
		b[i]=i*(i+alfa);
	}
	em[0]=FLT_MIN;
	em[2]=FLT_EPSILON;
	em[4]=6*n;
	allzerortpol(n,a,b,zer,em);
	free_real_vector(a,0);
	free_real_vector(b,0);
}

