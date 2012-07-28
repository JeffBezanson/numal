#include "../real.h"
void elmrow(int l, int u, int i, int j, real_t **a, real_t **b, real_t x)
{
	for (; l<=u; l++) a[i][l] += b[j][l]*x;
}
