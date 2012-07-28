#include "../real.h"


void allzerortpol(int n, real_t b[], real_t c[], real_t zer[],
						real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int qrivalsymtri(real_t [], real_t [], int, real_t []);
	void dupvec(int, int, int, real_t [], real_t []);
	int i;
	real_t nrm,*bb;

	bb=allocate_real_vector(1,n);
	nrm=fabs(b[0]);
	for (i=1; i<=n-2; i++)
		if (c[i]+fabs(b[i]) > nrm) nrm=c[i]+fabs(b[i]);
	if (n > 1)
		nrm = (nrm+1 >= c[n-1]+fabs(b[n-1])) ? nrm+1.0 :
					(c[n-1]+fabs(b[n-1]));
	em[1]=nrm;
	for (i=n; i>=1; i--) zer[i]=b[i-1];
	dupvec(1,n-1,0,bb,c);
	qrivalsymtri(zer,bb,n,em);
	free_real_vector(bb,1);
}
