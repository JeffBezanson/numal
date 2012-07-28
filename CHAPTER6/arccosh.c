#include "../real.h"


real_t arccosh(real_t x)
{
	return ((x <= 1.0) ? 0.0 : ((x > 1.0e10) ?
						0.69314718055995+log(x) :
						log(x+sqrt((x-1.0)*(x+1.0)))));
}
