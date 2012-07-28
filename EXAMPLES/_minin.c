#include <stdio.h>
#include <math.h>

real_t f(real_t x)
{
	int i;
	real_t s,temp;

	s=0.0;
	for (i=1; i<=20; i++) {
		temp=(i*2-5)/(x-i*i);
		s += temp*temp;
	}
	return s;
}

real_t tol(real_t x)
{
	return (fabs(x)*1.0e-6+1.0e-6);
}

void main ()
{
	real_t minin(real_t *, real_t *, real_t *, real_t (*)(real_t),
					real_t (*)(real_t));
	real_t m,x,a,b;

	a=1.0+tol(1.0);
	b=4.0-tol(4.0);
	m=minin(&x,&a,&b,f,tol);
	printf("Minimum is  %e\nFor x is  %e\nin the interval with "
			"endpoints  %e   %e",m,x,a,b);
}

