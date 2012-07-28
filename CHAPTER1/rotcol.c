#include "../real.h"
void rotcol(int l, int u, int i, int j, real_t **a, real_t c, real_t s)
{
	real_t x, y;

	for (; l<=u; l++) {
		x=a[l][i];
		y=a[l][j];
		a[l][i]=x*c+y*s;
		a[l][j]=y*c-x*s;
	}
}
