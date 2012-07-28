#include <stdio.h>
void main ()
{
	void bessk(real_t, int, real_t []);
	int j;
	real_t x,k[3];

	printf("BESSK delivers:\n");
	x=0.5;
	for (j=1; j<=4; j++) {
		bessk(x,2,k);
		printf(" %4.1f   %e   %e   %e\n",x,k[0],k[1],k[2]);
		x += 0.5;
	}
}

