#include "../real.h"
void dupcolvec(int l, int u, int j, real_t **a, real_t b[])
{
	for (; l<=u; l++) a[l][j]=b[l];
}

