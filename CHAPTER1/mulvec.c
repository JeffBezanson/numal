#include "../real.h"
void mulvec(int l, int u, int shift, real_t a[], real_t b[], real_t x)
{
	for (; l<=u; l++) a[l]=b[l+shift]*x;
}
