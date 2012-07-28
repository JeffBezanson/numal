#include "../real.h"


real_t logoneplusx(real_t x)
{
	real_t y,z;

	if (x == 0.0)
		return 0.0;
	else if (x < -0.2928 || x > 0.4142)
		return log(1.0+x);
	else {
		z=x/(x+2.0);
		y=z*z;
		return z*(2.0+y*(0.666666666663366+y*(0.400000001206045+y*
					(0.285714091590488+y*(0.22223823332791+y*
					(0.1811136267967+y*0.16948212488))))));
	}
}
