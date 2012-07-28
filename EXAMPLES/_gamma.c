#include <stdio.h>
void main ()
{
	real_t gamma(real_t);
	int i;
	static real_t x[4]={-8.5, 0.25, 1.5, 22.0};
	printf("GAMMA delivers:\n");
	for (i=0; i<=3; i++)
		printf(" %6.2f   %+e\n",x[i],gamma(x[i]));
}

