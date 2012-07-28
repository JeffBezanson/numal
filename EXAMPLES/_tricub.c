#include <math.h>
#include <stdio.h>

real_t e(real_t x, real_t y)
{
	return cos(x)*cos(y);
}

void main ()
{
	real_t tricub(real_t, real_t, real_t, real_t, real_t, real_t,
					real_t (*)(real_t,real_t), real_t, real_t);
	int i;
	real_t pi,acc;

	printf("TRICUB delivers:\n");
	pi=3.14159265359;
	acc=1.0;
	for (i=0; i<=5; i++) {
		acc *= 1.0e-1;
		printf(" %9.1e   %e\n",acc,
				tricub(0.0,0.0,0.0,pi/2.0,pi/2.0,pi/2.0,e,acc,acc));
	}
}

