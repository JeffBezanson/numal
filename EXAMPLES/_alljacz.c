#include <stdio.h>
void main ()
{
	void alljaczer(int, real_t, real_t, real_t []);
	real_t x[4];

	alljaczer(3,-0.5,-0.5,x);
	printf("Delivers:\n  %13.6e %13.6e %13.6e\n",x[1],x[2],x[3]);
}

