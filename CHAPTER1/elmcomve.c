#include "../real.h"
void elmcomveccol(int l, int u, int j, real_t ar[], real_t ai[],
						real_t **br, real_t **bi, real_t xr, real_t xi)
{
	void elmveccol(int, int, int, real_t [], real_t **, real_t);

	elmveccol(l,u,j,ar,br,xr);
	elmveccol(l,u,j,ar,bi,-xi);
	elmveccol(l,u,j,ai,br,xi);
	elmveccol(l,u,j,ai,bi,xr);
}
