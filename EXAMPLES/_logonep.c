#include <stdio.h>
void main ()
{
	real_t logoneplusx(real_t);
	int i;
	real_t x;

	printf("LOGONEPLUSX delivers:\n");
	x=1.0e-1;
	for (i=1; i<=2; i++) {
		printf("  %e\n",logoneplusx(2.0*exp(x/2.0)*sinh(x/2.0)));
		x *= 1.0e-9;
	}
}

