#include <math.h>
#include <stdio.h>

real_t a(real_t x)
{
	return (sin(x));
}

void main ()
{
	real_t qadrat(real_t *, real_t, real_t, real_t (*)(real_t), real_t []);
	real_t t,q,e[4];

	e[1]=e[2]=1.0e-6;
	q=qadrat(&t,0.0,3.141592653589,a,e);
	printf("Delivers:  %13.6e  %3.0f\n",q,e[3]);
}

