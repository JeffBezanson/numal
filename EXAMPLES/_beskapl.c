#include <stdio.h>
void main ()
{
	void besskaplusn(real_t, real_t, int, real_t []);
	real_t kan[3];

	besskaplusn(0.0,1.0,2,kan);
	printf("BESSKAPLUSN delivers:\n %e  %e  %e\n",
			kan[0],kan[1],kan[2]);
}

