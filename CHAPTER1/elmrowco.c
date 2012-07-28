#include "../real.h"
void elmrowcol(int l, int u, int i, int j, real_t **a, real_t **b, real_t x)
{
	for (; l<=u; l++) a[i][l] += b[l][j]*x;
}
