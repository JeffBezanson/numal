#include "../real.h"
void gsswtssym(int n, real_t zer[], real_t c[], real_t w[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void allortpolsym(int, real_t, real_t [], real_t []);
	int i,twoi,low,up;
	real_t s,*p;

	p=allocate_real_vector(0,n-1);
	low=1;
	up=n;
	while (low < up) {
		allortpolsym(n-1,zer[low],c,p);
		s=p[n-1]*p[n-1];
		for (i=n-1; i>=1; i--) s=s/c[i]+p[i-1]*p[i-1];
		w[low]=1.0/s;
		low++;
		up--;
	}
	if (low == up) {
		s=1.0;
		for (twoi=n-1; twoi>=2; twoi -= 2) s=s*c[twoi-1]/c[twoi]+1.0;
		w[low]=1.0/s;
	}
	free_real_vector(p,0);
}
