#include <stdio.h>
void main ()
{
	void eialpha(real_t, int, real_t []);
	int k;
	real_t a[6];

	printf("EIALPHA delivers:\n");
	eialpha(0.25,5,a);
	for (k=0; k<=5; k++) {
		printf(" %2d   %e\n",k,a[k]);
	}
}

