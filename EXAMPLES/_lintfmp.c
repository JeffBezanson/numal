#include <stdio.h>
void main ()
{
	void lintfmpol(real_t, real_t, int, real_t []);
	real_t a[3] = {1.0, 2.0, 3.0};

	printf("              a[0]   a[1]   a[2]\ninput:"
			"     %7.2f%7.2f%7.2f",a[0],a[1],a[2]);
	lintfmpol(2.0,3.0,2,a);
	printf("\nlintfmpol: %7.2f%7.2f%7.2f  (power sum in y)",
			a[0],a[1],a[2]);
	lintfmpol(0.5,-1.5,2,a);
	printf("\nlintfmpol: %7.2f%7.2f%7.2f  (power sum in x)",
			a[0],a[1],a[2]);
}

