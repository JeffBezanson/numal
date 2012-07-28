#include "../real.h"
real_t infnrmmat(int lr, int ur, int lc, int uc, int *kr, real_t **a)
{
	real_t onenrmrow(int, int, int, real_t **);
	real_t r, max;

	max=0.0;
	*kr=lr;
	for (; lr<=ur; lr++) {
		r=onenrmrow(lc,uc,lr,a);
		if (r > max) {
			max=r;
			*kr=lr;
		}
	}
	return (max);
}
