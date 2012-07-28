#include "../real.h"


real_t bessi1(real_t x)
{
	if (x == 0.0) return 0.0;
	if (fabs(x) <= 15.0) {
		real_t z,denominator,numerator;
		z=x*x;
		denominator=z*(z-0.222583674000860e4)+0.136293593052499e7;
		numerator=(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*
			0.207175767232792e-26+0.257091905584414e-23)+
			0.306279283656135e-20)+0.261372772158124e-17)+
			0.178469361410091e-14)+0.963628891518450e-12)+
			0.410068906847159e-9)+0.135455228841096e-6)+
			0.339472890308516e-4)+0.624726195127003e-2)+
			0.806144878821295e0)+0.682100567980207e2)+
			0.341069752284422e4)+0.840705772877836e5)+
			0.681467965262502e6);
		return x*(numerator/denominator);
	} else {
		real_t nonexpbessi1(real_t);
		return exp(fabs(x))*nonexpbessi1(x);
	}
}
