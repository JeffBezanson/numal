#include "../real.h"
void fulsymmatvec(int lr, int ur, int lc, int uc,
						real_t a[], real_t b[], real_t c[])
{
	real_t symmatvec(int, int, int, real_t [], real_t []);

	for (; lr<=ur; lr++) c[lr]=symmatvec(lc,uc,lr,a,b);
}
