#include <stdio.h>
void main ()
{
	real_t arcsinh(real_t);
	printf("ARCSINH delivers:\n");
	printf("  %e\n",arcsinh(sinh(0.01)));
	printf("  %e\n",arcsinh(sinh(0.05)));
	printf("  %e\n",sinh(arcsinh(0.05)));
	printf("  %e\n",sinh(arcsinh(0.01)));
}

