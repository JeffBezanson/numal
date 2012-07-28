#include "../real.h"
void elmvecrow(int l, int u, int i, real_t a[], real_t **b, real_t x)
{
	for (; l<=u; l++) a[l] += b[i][l]*x;
}
