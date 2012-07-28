#include <stdio.h>
void main ()
{
	real_t arctanh(real_t);
	printf("ARCTANH delivers:\n");
	printf("  %e\n",arctanh(tanh(0.01)));
	printf("  %e\n",arctanh(tanh(0.05)));
	printf("  %e\n",tanh(arctanh(0.05)));
	printf("  %e\n",tanh(arctanh(0.01)));
}

