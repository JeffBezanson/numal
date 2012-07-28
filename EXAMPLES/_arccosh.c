#include <stdio.h>
void main ()
{
	real_t arccosh(real_t);
	printf("ARCCOSH delivers:\n");
	printf("  %e\n",arccosh(cosh(0.01)));
	printf("  %e\n",arccosh(cosh(0.05)));
	printf("  %e\n",cosh(arccosh(1.01)));
	printf("  %e\n",cosh(arccosh(1.05)));
}

