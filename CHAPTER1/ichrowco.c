#include "../real.h"
void ichrowcol(int l, int u, int i, int j, real_t **a)
{
	real_t r;

	for (; l<=u; l++) {
		r=a[i][l];
		a[i][l]=a[l][j];
		a[l][j]=r;
	}
}
