#include <math.h>
#include <stdio.h>
void main ()
{
	void spherbessi(real_t, int, real_t []);
	void nonexpspherbessi(real_t, int, real_t []);
	int n;
	real_t x,expx,i1[4],i2[4];

	printf("SPHERBESSI and NONEXPSPHERBESSI deliver:\n");
	x=1.0;
	expx=exp(x);
	n=3;
	spherbessi(x,n,i1);
	nonexpspherbessi(x,n,i2);
	for (n=0; n<=3; n++)
		printf("  %d   %e  %e\n",n,i1[n],i2[n]*expx);
}

