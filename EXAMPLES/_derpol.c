#include <stdio.h>
void main ()
{
	void derpol(int, int, real_t, real_t []);
	real_t a[4] = {-1.0, 1.0, -2.0, 3.0};

	derpol(3,3,1.0,a);
	printf("The 0-th until and including the 3rd derivatives:"
			"  %-4.0f%-4.0f%-4.0f%-4.0f",a[0],a[1],a[2],a[3]);
}

