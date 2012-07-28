#include <stdio.h>
void main ()
{
	void bessya01(real_t, real_t, real_t *, real_t *);
	real_t p,q;

	bessya01(0.0,1.0,&p,&q);
	printf("BESSYA01 delivers:  %e   %e\n",p,q);
}

