#include <stdio.h>
void main ()
{
	void ibpplusn(real_t, real_t, real_t, int, real_t, real_t []);
	real_t isubx[3];

	ibpplusn(0.3,0.4,1.5,2,1.0e-6,isubx);
	printf("IBPPLUSN delivers:\n %e   %e   %e\n",
			isubx[0],isubx[1],isubx[2]);
}

