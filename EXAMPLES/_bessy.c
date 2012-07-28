#include <stdio.h>
void main ()
{
	void bessy(real_t, int, real_t []);
	real_t y[3];

	bessy(1.0,2,y);
	printf("BESSY delivers:\n  %e   %e   %e\n",y[0],y[1],y[2]);
}

