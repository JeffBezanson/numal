#include "../real.h"
#include <stdlib.h>


void inisymd(int lr, int ur, int shift, real_t a[], real_t x)
{
	shift=abs(shift);
	ur += shift+1;
	shift += lr;
	lr += ((shift-3)*shift)/2;
	lr += shift;
	while (shift < ur) {
		a[lr]=x;
		shift++;
		lr += shift;
	}
}
