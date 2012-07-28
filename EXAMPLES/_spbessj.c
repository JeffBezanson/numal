#include <stdio.h>
void main ()
{
	void spherbessj(real_t, int, real_t []);
	int n;
	real_t x,j[3];

	x=1.5;
	n=2;
	spherbessj(x,n,j);
	printf("SPHERBESSJ delivers:\n X = %3.1f    N = %d\n"
			" %e  %e  %e\n",x,n,j[0],j[1],j[2]);
}

