#include <math.h>
#include <stdio.h>

real_t b(real_t x, real_t y)
{
	return x+y;
}

real_t fxy(real_t x, real_t y)
{
	return 1.0-2.0*(x*x+y);
}

void main ()
{
	void rk4a(real_t *, real_t, real_t (*)(real_t, real_t),
				real_t *, real_t, real_t (*)(real_t, real_t),
				real_t [], real_t [], int, int, int);
	real_t x,y,d[5],e[6];

	e[0]=e[1]=e[2]=e[3]=e[4]=e[5]=1.0e-4;
	rk4a(&x,0.0,b,&y,0.0,fxy,e,d,1,1,1);
	printf("x =  %e        Exactly : 2.00000\ny = %e\n"
				"y-x*(1-x) = %e\n",x,y,y-x*(1-x));
}

