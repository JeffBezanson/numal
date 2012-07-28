#include "../real.h"
void fulmatvec(int lr, int ur, int lc, int uc,
					real_t **a, real_t b[], real_t c[])
{
	real_t matvec(int, int, int, real_t **, real_t []);

	for (; lr<=ur; lr++) c[lr]=matvec(lc,uc,lr,a,b);
}
