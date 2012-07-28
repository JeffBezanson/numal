#include <stdio.h>
void main ()
{
	real_t p,q;
	void incomgam(real_t, real_t, real_t *, real_t *, real_t, real_t);
	printf("INCOMGAM delivers:\n");
	incomgam(3.0,4.0,&p,&q,1.0*2.0*3.0,1.0e-6);
	printf("  KLGAM and GRGAM are : %e  %e\n",p,q);
}

