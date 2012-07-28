#include "../real.h"
void ichcol(int l, int u, int i, int j, real_t **a)
{
	real_t r;

	for (; l<=u; l++) {
		r=a[l][i];
		a[l][i]=a[l][j];
		a[l][j]=r;
	}
}
