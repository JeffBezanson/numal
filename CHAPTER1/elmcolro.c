#include "../real.h"
void elmcolrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x)
{
	for (; l<=u; l++) a[l][i] += b[j][l]*x;
}
