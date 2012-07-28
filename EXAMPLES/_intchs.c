#include <math.h>
#include <stdio.h>
void main ()
{
	void intchs(int, real_t [], real_t []);
	real_t b[5], a[4] = {1.0, 0.5, 0.2, 0.1};

	intchs(3,a,b);
	printf("INTCHS delivers:%8.4f%8.4f%8.4f%8.4f",b[1],b[2],b[3],b[4]);
}

