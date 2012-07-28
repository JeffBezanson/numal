#include "../real.h"


void spherbessk(real_t x, int n, real_t k[])
{
	void nonexpspherbessk(real_t, int, real_t []);
	real_t expx;
	expx=exp(-x);
	nonexpspherbessk(x,n,k);
	for (; n>=0; n--) k[n] *= expx;
}
