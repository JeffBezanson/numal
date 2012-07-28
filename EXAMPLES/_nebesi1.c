#include <stdio.h>
void main ()
{
	real_t nonexpbessi1(real_t);
	real_t x;

	x=1.0;
	printf("NONEXPBESSI1 delivers:  %2.0f    %e\n",
			x,nonexpbessi1(x));
}

