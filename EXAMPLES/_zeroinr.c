#include <math.h>
#include <stdio.h>

real_t fx(real_t  x)
{
	return exp(-x*3.0)*(x-1.0)+x*x*x;
}

real_t tolx(real_t  x)
{
	return fabs(x)*1.0e-6+1.0e-6;
}

void main ()
{
	int zeroinrat(real_t *, real_t *, real_t (*)(real_t), real_t (*)(real_t));
	real_t x,y;

	x=0.0;
	y=1.0;
	if (zeroinrat(&x,&y,fx,tolx))
		printf("Calculated zero: %e\n",x);
	else
		printf("No zero found.");
}

