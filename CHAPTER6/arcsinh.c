#include "../real.h"


real_t arcsinh(real_t x)
{
	real_t logoneplusx(real_t);
	real_t y;

	if (fabs(x) > 1.0e10)
		return ((x > 0.0) ? 0.69314718055995+log(fabs(x)) :
									-0.69314718055995+log(fabs(x)));
	else {
		y=x*x;
		return ((x == 0.0) ? 0.0 :	((x > 0.0) ?
					logoneplusx(fabs(x)+y/(1.0+sqrt(1.0+y))) :
					-logoneplusx(fabs(x)+y/(1.0+sqrt(1.0+y)))));
	}
}
