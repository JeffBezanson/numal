#include <stdio.h>

real_t a(real_t x)
{
	return sqrt(1.0-x*x);
}

void main ()
{
	void reccof(int, int, real_t *, real_t (*)(real_t), real_t [],
					real_t [], real_t [], int);
	real_t x,b[3],c[3],l[3];

	reccof(2,200,&x,a,b,c,l,1);
	printf("Delivers: %7.3f %7.3f\n",c[1],c[2]);
}

