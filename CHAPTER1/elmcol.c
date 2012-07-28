#include "../real.h"
void elmcol(int l, int u, int i, int j, real_t **a, real_t **b, real_t x)
{
	for (; l<=u; l++) a[l][i] += b[l][j]*x;
}
