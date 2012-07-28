#include "../real.h"


real_t comeucnrm(real_t **ar, real_t **ai, int lw, int n)
{
	real_t mattam(int, int, int, int, real_t **, real_t **);
	int i,l;
	real_t r;

	r=0.0;
	for (i=1; i<=n; i++) {
		l=(i>lw) ? i-lw : 1;
		r += mattam(l,n,i,i,ar,ar)+mattam(l,n,i,i,ai,ai);
	}
	return (sqrt(r));
}
