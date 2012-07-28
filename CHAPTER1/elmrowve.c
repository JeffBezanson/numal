#include "../real.h"
void elmrowvec(int l, int u, int i, real_t **a, real_t b[], real_t x)
{
	for (; l<=u; l++) a[i][l] += b[l]*x;
}
