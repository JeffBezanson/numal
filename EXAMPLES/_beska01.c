#include <stdio.h>
void main ()
{
	void besska01(real_t, real_t, real_t *, real_t *);
	real_t p,q;

	besska01(0.0,1.0,&p,&q);
	printf("BESSKA01 delivers:  %e  %e\n",p,q);
}

