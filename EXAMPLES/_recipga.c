#include <stdio.h>
void main ()
{
	real_t recipgamma(real_t, real_t *, real_t *);
	real_t x,odd,even;
	printf("RECIPGAMMA delivers:\n");
	x=recipgamma(0.4,&odd,&even);
	printf(" 0.4   %e   %e   %e\n",x,odd,even);
	x=recipgamma(0.0,&odd,&even);
	printf(" 0.0   %e   %e   %e\n",x,odd,even);
}

