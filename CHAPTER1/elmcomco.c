#include "../real.h"
void elmcomcol(int l, int u, int i, int j, real_t **ar, real_t **ai,
					real_t **br, real_t **bi, real_t xr, real_t xi)
{
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);

	elmcol(l,u,i,j,ar,br,xr);
	elmcol(l,u,i,j,ar,bi,-xi);
	elmcol(l,u,i,j,ai,br,xi);
	elmcol(l,u,i,j,ai,bi,xr);
}
