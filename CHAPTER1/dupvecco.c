#include "../real.h"
void dupveccol(int l, int u, int j, real_t a[], real_t **b)
{
	for (; l<=u; l++) a[l]=b[l][j];
}

