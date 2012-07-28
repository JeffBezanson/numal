#include <stdio.h>
void main ()
{
	void nonexpbessi(real_t, int, real_t []);
	int j;
	real_t x,i[3];

	printf("NONEXPBESSI delivers:\n");
	x=0.5;
	for (j=1; j<=5; j++) {
		nonexpbessi(x,2,i);
		printf(" %4.1f   %e   %e   %e\n",x,i[0],i[1],i[2]);
		x += 0.5;
	}
}

