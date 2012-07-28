#include "../real.h"
void dupvec(int l, int u, int shift, real_t a[], real_t b[])
{
	for (; l<=u; l++) a[l]=b[l+shift];
}
