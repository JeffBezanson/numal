#include "../real.h"


real_t bessi0(real_t x)
{
	if (x == 0.0) return 1.0;
	if (fabs(x) <= 15.0) {
		real_t z,denominator,numerator;
		z=x*x;
		numerator=(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*
			0.210580722890567e-22+0.380715242345326e-19)+
			0.479440257548300e-16)+0.435125971262668e-13)+
			0.300931127112960e-10)+0.160224679395361e-7)+
			0.654858370096785e-5)+0.202591084143397e-2)+
			0.463076284721000e0)+0.754337328948189e2)+
			0.830792541809429e4)+0.571661130563785e6)+
			0.216415572361227e8)+0.356644482244025e9)+
			0.144048298227235e10);
		denominator=(z*(z*(z-0.307646912682801e4)+
			0.347626332405882e7)-0.144048298227235e10);
		return -numerator/denominator;
	} else {
		real_t nonexpbessi0(real_t);
		return exp(fabs(x))*nonexpbessi0(x);
	}
}
