#include "../real.h"
void rotrow(int l, int u, int i, int j, real_t **a, real_t c, real_t s)
{
	real_t x, y;

	for (; l<=u; l++) {
		x=a[i][l];
		y=a[j][l];
		a[i][l]=x*c+y*s;
		a[j][l]=y*c-x*s;
	}
}
