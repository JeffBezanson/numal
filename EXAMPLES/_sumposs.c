#include <stdio.h>

real_t ai(real_t i)
{
	return 1.0/(i*i);
}

void main ()
{
	real_t sumposseries(real_t (*)(real_t), int, real_t, int, int, int);

	printf("SUMPOSSERIES delivers:  %e\n",
			sumposseries(ai,100,1.0e-6,8,100,10));
}

