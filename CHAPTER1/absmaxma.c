#include "../real.h"
real_t absmaxmat(int lr, int ur, int lc, int uc, int *i, int *j, real_t **a)
{
	real_t infnrmcol(int, int, int, int *, real_t **);
	int ii;
	real_t r, max;

	max=0.0;
	*i=lr;
	*j=lc;
	for (; lc<=uc; lc++) {
		r=infnrmcol(lr,ur,lc,&ii,a);
		if (r > max) {
			max=r;
			*i=ii;
			*j=lc;
		}
	}
	return (max);
}
