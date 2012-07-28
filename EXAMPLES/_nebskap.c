#include <stdio.h>
void main ()
{
	void nonexpbesskaplusn(real_t, real_t, int, real_t []);
	real_t kan[3];

	nonexpbesskaplusn(0.0,5.0,2,kan);
	printf("NONEXPBESSKAPLUSN delivers:\n"
			" %e  %e  %e\n",kan[0],kan[1],kan[2]);
}

