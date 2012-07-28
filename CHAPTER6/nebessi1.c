#include "../real.h"


real_t nonexpbessi1(real_t x)
{
	if (x == 0.0) return 0.0;
	if (fabs(x) > 15.0) {
		int i,signx;
		real_t br,br1,br2,z,z2,sqrtx,numerator,denominator;
		static real_t ar1[4]={0.1494052814740e1, -0.362026420242263e3,
			0.220549722260336e5, -0.408928084944275e6};
		static real_t ar2[4]={1.0, -0.631003200551590e3,
			0.496811949533398e5, -0.100425428133695e7};
		signx = (x > 0.0) ? 1 : -1;
		x=fabs(x);
		sqrtx=sqrt(x);
		z=30.0/x-1.0;
		z2=z+z;
		br1=br2=0.0;
		for (i=0; i<=3; i++) {
			br=z2*br1-br2+ar1[i];
			br2=br1;
			br1=br;
		}
		numerator=z*br1-br2+0.102776692371524e7;
		br1=br2=0.0;
		for (i=0; i<=3; i++) {
			br=z2*br1-br2+ar2[i];
			br2=br1;
			br1=br;
		}
		denominator=z*br1-br2+0.26028876789105e7;
		return ((numerator/denominator)/sqrtx)*signx;
	} else {
		real_t bessi1(real_t);
		return exp(-fabs(x))*bessi1(x);
	}
}
