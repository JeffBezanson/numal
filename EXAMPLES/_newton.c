#include <stdio.h>
void main ()
{
	void newton(int, real_t [], real_t []);
	real_t x[3],f[3];

	x[0]=0.0;  x[1]=0.5;  x[2]=1.0;
	f[0]=1.0;  f[1]=f[2]=0.0;
	newton(2,x,f);
	printf("The Newton coefficients are:\n %e   %e   %e",
			f[0],f[1],f[2]);
}

