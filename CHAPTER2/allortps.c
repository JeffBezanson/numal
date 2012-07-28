#include "../real.h"
void allortpolsym(int n, real_t x, real_t c[], real_t p[])
{
	int k;
	real_t r,s,h;

	if (n == 0) {
		p[0]=1.0;
		return;
	}
	r=p[1]=x;
	s=p[0]=1.0;
	for (k=2; k<=n; k++) {
		h=r;
		p[k]=r=x*r-c[k-1]*s;
		s=h;
	}
}
