#include "../real.h"
void comrowcst(int l, int u, int i,
					real_t **ar, real_t **ai, real_t xr, real_t xi)
{
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);

	for (; l<=u; l++)
		commul(ar[i][l],ai[i][l],xr,xi,&ar[i][l],&ai[i][l]);
}
