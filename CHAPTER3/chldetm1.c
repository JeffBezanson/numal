#include "../real.h"
real_t chldeterm1(real_t a[], int n)
{
	int k,kk;
	real_t d;

	d=1.0;
	kk=0;
	for (k=1; k<=n; k++) {
		kk += k;
		d *= a[kk];
	}
	return (d*d);
}
