#include <stdio.h>
void main ()
{
	void ibqplusn(real_t, real_t, real_t, int, real_t, real_t []);
	real_t isubx[3];

	ibqplusn(0.3,1.4,0.5,2,1.0e-6,isubx);
	printf("IBQPLUSN delivers:\n %e   %e   %e\n",
			isubx[0],isubx[1],isubx[2]);
}

