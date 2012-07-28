#include "../real.h"
void dupvecrow(int l, int u, int i, real_t a[], real_t **b)
{
	for (; l<=u; l++) a[l]=b[i][l];
}
