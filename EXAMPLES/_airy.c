#include <math.h>
#include <stdio.h>
void main ()
{
	void airy(real_t, real_t *, real_t *, real_t *, real_t *,
				real_t *, int);
	real_t a,b,c,d,e;

	airy(9.654894,&a,&b,&c,&d,&e,1);
	printf("AIRY delivers:\n"
			" AI (9.654894) = %+e\n AID(9.654894) = %+e\n"
			" BI (9.654894) = %+e\n BID(9.654894) = %+e\n",
			a*exp(-e),b*exp(-e),c*exp(e),d*exp(e));
}

