#include "../real.h"


void errorfunction(real_t x, real_t *erf, real_t *erfc)
{
	if (x > 26.0) {
		*erf = 1.0;
		*erfc = 0.0;
		return;
	} else if (x < -5.5) {
		*erf = -1.0;
		*erfc = 2.0;
		return;
	} else {
		real_t nonexperfc(real_t);
		real_t absx,c,p,q;
		absx=fabs(x);
		if (absx <= 0.5) {
			c=x*x;
			p=((-0.356098437018154e-1*c+0.699638348861914e1)*c+
					0.219792616182942e2)*c+0.242667955230532e3;
			q=((c+0.150827976304078e2)*c+0.911649054045149e2)*c+
					0.215058875869861e3;
			*erf = x*p/q;
			*erfc = 1.0-(*erf);
		} else {
			*erfc = exp(-x*x)*nonexperfc(absx);
			*erf = 1.0-(*erfc);
			if (x < 0.0) {
				*erf = -(*erf);
				*erfc = 2.0-(*erfc);
			}
		}
	}
}
