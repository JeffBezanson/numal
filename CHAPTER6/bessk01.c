#include "../real.h"


void bessk01(real_t x, real_t *k0, real_t *k1)
{
	if (x <= 1.5) {
		int k;
		real_t c,d,r,sum0,sum1,t,term,t0,t1;
		sum0=d=log(2.0/x)-0.5772156649015328606;
		sum1 = c = -1.0-2.0*d;
		r=term=1.0;
		t=x*x/4.0;
		k=1;
		do {
			term *= t*r*r;
			d += r;
			c -= r;
			r=1.0/(k+1);
			c -= r;
			t0=term*d;
			t1=term*c*r;
			sum0 += t0;
			sum1 += t1;
			k++;
		} while (fabs(t0/sum0)+fabs(t1/sum1) > 1.0e-15);
		*k0 = sum0;
		*k1 = (1.0+t*sum1)/x;
	} else {
		void nonexpbessk01(real_t, real_t *, real_t *);
		real_t expx;
		expx=exp(-x);
		nonexpbessk01(x,k0,k1);
		*k1 *= expx;
		*k0 *= expx;
	}
}
