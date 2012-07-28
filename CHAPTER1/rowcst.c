#include "../real.h"
void rowcst(int l, int u, int i, real_t **a, real_t x)
{
	for (; l<=u; l++) a[i][l] *= x;
}
