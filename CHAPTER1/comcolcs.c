#include "../real.h"
void comcolcst(int l, int u, int j,
					real_t **ar, real_t **ai, real_t xr, real_t xi)
{
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);

	for (; l<=u; l++)
		commul(ar[l][j],ai[l][j],xr,xi,&ar[l][j],&ai[l][j]);
}
