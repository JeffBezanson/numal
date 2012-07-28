#include <math.h>
#include <stdio.h>

real_t f(real_t  x)
{
	return exp(-x*3.0)*(x-1.0)+x*x*x;
}

real_t df(real_t  x)
{
	return exp(-x*3.0)*(-3.0*x+4.0)+3.0*x*x;
}

real_t tolx(real_t  x)
{
	return fabs(x)*1.0e-6+1.0e-6;
}

void main ()
{
	int zeroinder(real_t *, real_t *, real_t (*)(real_t),
						real_t (*)(real_t), real_t (*)(real_t));
	real_t x,y;

	x=0.0;
	y=1.0;
	if (zeroinder(&x,&y,f,df,tolx))
		printf("Calculated zero and function value:\n  %e     %e"
				"\nOther straddling approximation and function value:"
				"\n  %e    %e\n",x,f(x),y,f(y));
	else
		printf("No zero found.");
}

