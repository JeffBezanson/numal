#include "../real.h"
void colcst(int l, int u, int j, real_t **a, real_t x)
{
	for (; l<=u; l++) a[l][j] *= x;
}
