#include "../real.h"
void duprowvec(int l, int u, int i, real_t **a, real_t b[])
{
	for (; l<=u; l++) a[i][l]=b[l];
}
