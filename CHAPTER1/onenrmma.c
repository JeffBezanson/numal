#include "../real.h"
real_t onenrmmat(int lr, int ur, int lc, int uc, int *kc, real_t **a)
{
	real_t onenrmcol(int l, int u, int j, real_t **a);
	real_t r, max;

	max=0.0;
	*kc=lc;
	for (; lc<=uc; lc++) {
		r=onenrmcol(lr,ur,lc,a);
		if (r > max) {
			max=r;
			*kc=lc;
		}
	}
	return (max);
}
