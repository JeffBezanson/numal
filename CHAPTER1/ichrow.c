#include "../real.h"
void ichrow(int l, int u, int i, int j, real_t **a)
{
	real_t r;

	for (; l<=u; l++) {
		r=a[i][l];
		a[i][l]=a[j][l];
		a[j][l]=r;
	}
}
