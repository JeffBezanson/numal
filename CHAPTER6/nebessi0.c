#include "../real.h"


real_t nonexpbessi0(real_t x)
{
	if (x == 0.0) return 1.0;
	if (fabs(x) <= 15.0) {
		real_t bessi0(real_t);
		return exp(-fabs(x))*bessi0(x);
	} else {
		int i;
		real_t sqrtx,br,br1,br2,z,z2,numerator,denominator;
		static real_t ar1[4]={0.2439260769778, -0.115591978104435e3,
			0.784034249005088e4, -0.143464631313583e6};
		static real_t ar2[4]={1.0, -0.325197333369824e3,
			0.203128436100794e5, -0.361847779219653e6};
		x=fabs(x);
		sqrtx=sqrt(x);
		br1=br2=0.0;
		z=30.0/x-1.0;
		z2=z+z;
		for (i=0; i<=3; i++) {
			br=z2*br1-br2+ar1[i];
			br2=br1;
			br1=br;
		}
		numerator=z*br1-br2+0.346519833357379e6;
		br1=br2=0.0;
		for (i=0; i<=3; i++) {
			br=z2*br1-br2+ar2[i];
			br2=br1;
			br1=br;
		}
		denominator=z*br1-br2+0.865665274832055e6;
		return (numerator/denominator)/sqrtx;
	}
}
