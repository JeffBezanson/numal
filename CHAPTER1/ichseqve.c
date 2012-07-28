#include "../real.h"
void ichseqvec(int l, int u, int il, int shift, real_t a[])
{
	real_t r;

	for (; l<=u; l++) {
		r=a[il];
		a[il]=a[l+shift];
		a[l+shift]=r;
		il += l;
	}
}
