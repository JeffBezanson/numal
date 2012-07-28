#include "../real.h"
void dupmat(int l, int u, int i, int j, real_t **a, real_t **b)
{
	int k;

	for (; l<=u; l++)
		for (k=i; k<=j; k++) a[l][k]=b[l][k];
}
