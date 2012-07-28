#include "../real.h"


void bessk(real_t x, int n, real_t k[])
{
	void bessk01(real_t, real_t *, real_t *);
	int i;
	real_t k0,k1,k2;

	bessk01(x,&k0,&k1);
	k[0]=k0;
	if (n > 0) k[1]=k1;
	x=2.0/x;
	for (i=2; i<=n; i++) {
		k[i]=k2=k0+x*(i-1)*k1;
		k0=k1;
		k1=k2;
	}
}
