#include <stdio.h>
void main ()
{
	void errorfunction(real_t, real_t *, real_t *);
	real_t nonexperfc(real_t);
	real_t erf,erfc,p;

	errorfunction(1.0,&erf,&erfc);
	p=nonexperfc(100.0);
	printf("ERRORFUNCTION and NONEXPERFC deliver:\n\n"
			"ERF(1)  = %e\nERFC(1) = %e\nNONEXPERFC(100) = %e\n",
			erf,erfc,p);
}

