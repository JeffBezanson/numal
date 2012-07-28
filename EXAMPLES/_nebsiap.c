#include <stdio.h>
void main ()
{
	void nonexpbessiaplusn(real_t, real_t, int, real_t []);
	real_t a,x,ia[3];

	a=0.25;
	x=2.0;
	nonexpbessiaplusn(a,x,2,ia);
	printf("NONEXPBESSIAPLUSN delivers:\n A = %4.2f    X = %4.2f\n"
			" %e   %e   %e\n",a,x,ia[0],ia[1],ia[2]);
}

