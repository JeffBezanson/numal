#include <stdio.h>

real_t a(int i)
{
	return (pow(-1,i)/((i+1)*(i+1)));
}

void main ()
{
	real_t euler(real_t (*)(int), real_t, int);

	printf("Delivers:  %13.6e\n",euler(a,1.0e-6,100));
}

