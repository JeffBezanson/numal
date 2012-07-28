#include "../real.h"
void ichseq(int l, int u, int il, int shift, real_t a[])
{
	real_t r;

	for (; l<=u; l++) {
		r=a[il];
		a[il]=a[il+shift];
		a[il+shift]=r;
		il += l;
	}
}
