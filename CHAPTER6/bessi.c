#include "../real.h"


void bessi(real_t x, int n, real_t i[])
{
	if (x == 0.0) {
		i[0]=1.0;
		for (; n>=1; n--) i[n]=0.0;
	} else {
		void nonexpbessi(real_t, int, real_t []);
		real_t expx;
		expx=exp(fabs(x));
		nonexpbessi(x,n,i);
		for (; n>=0; n--) i[n] *= expx;
	}
}
