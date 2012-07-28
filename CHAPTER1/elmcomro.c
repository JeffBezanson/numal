#include "../real.h"
void elmcomrowvec(int l, int u, int i, real_t **ar, real_t **ai,
						real_t br[], real_t bi[], real_t xr, real_t xi)
{
	void elmrowvec(int, int, int, real_t **, real_t [], real_t);

	elmrowvec(l,u,i,ar,br,xr);
	elmrowvec(l,u,i,ar,bi,-xi);
	elmrowvec(l,u,i,ai,br,xi);
	elmrowvec(l,u,i,ai,bi,xr);
}
