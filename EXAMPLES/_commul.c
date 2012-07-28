#include <stdio.h>
void main ()
{
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);
	real_t r,i;

	commul(0.1,0.2,0.3,0.4,&r,&i);
	printf("(.1+.2i)*(.3+.4*i) = %-5.2f+%-4.2f*i",r,i);
}

