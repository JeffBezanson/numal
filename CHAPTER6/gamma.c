#include "../real.h"


real_t numal_gamma(real_t x)
{
	real_t recipgamma(real_t, real_t *, real_t *);
	real_t loggamma(real_t);
	int inv;
	real_t y,s,f,g,odd,even;

	if (x < 0.5) {
		y=x-floor(x/2.0)*2;
		s = 3.14159265358979;
		if (y >= 1.0) {
			s = -s;
			y=2.0-y;
		}
		if (y >= 0.5) y=1.0-y;
		inv=1;
		x=1.0-x;
		f=s/sin(3.14159265358979*y);
	} else
		inv=0;
	if (x > 22.0)
		g=exp(loggamma(x));
	else {
		s=1.0;
		while (x > 1.5) {
			x=x-1.0;
			s *= x;
		}
		g=s/recipgamma(1.0-x,&odd,&even);
	}
	return (inv ? f/g : g);
}
