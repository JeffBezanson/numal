#include <stdio.h>
void main ()
{
	void inverseerrorfunction(real_t, real_t, real_t *);
	real_t inverf1,inverf2;

	inverseerrorfunction(0.6,0.0,&inverf1);
	inverseerrorfunction(1.0,1.0e-30,&inverf2);
	printf("INVERSEERRORFUNCTION delivers:\n\n"
			" X = 0.6,  INVERF = %e\n X = 1.0,  INVERF = %e\n",
			inverf1,inverf2);
}

