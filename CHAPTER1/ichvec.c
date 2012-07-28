#include "../real.h"
void ichvec(int l, int u, int shift, real_t a[])
{
	real_t r;

	for (; l<=u; l++) {
		r=a[l];
		a[l]=a[l+shift];
		a[l+shift]=r;
	}
}
