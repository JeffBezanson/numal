#include <math.h>
#include <stdio.h>

real_t fxy(real_t x, real_t y)
{
	return -y;
}

void main ()
{
	void rk1(real_t *, real_t, real_t, real_t *, real_t,
				real_t (*)(real_t, real_t), real_t [], real_t [], int);
	int first;
	real_t x,y,d[5],e[3];

	e[1]=e[2]=1.0e-4;
	first=1;
	rk1(&x,0.0,1.0,&y,1.0,fxy,e,d,first);
	printf("RK1 delivers:\n  x = %e\n  y = %e      yexact =  %e",
				x,y,exp(-x));
}

