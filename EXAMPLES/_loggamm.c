#include <stdio.h>
void main ()
{
	real_t loggamma(real_t);
	int i;
	static real_t x[5]={0.25, 1.5, 12.0, 15.0, 80.0};
	printf("LOGGAMMA delivers:\n");
	for (i=0; i<=4; i++)
		printf(" %6.2f   %+e\n",x[i],loggamma(x[i]));
}

