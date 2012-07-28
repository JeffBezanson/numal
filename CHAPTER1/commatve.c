#include "../real.h"
void commatvec(int l, int u, int i, real_t **ar, real_t **ai,
					real_t br[], real_t bi[], real_t *rr, real_t *ri)
{
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t mv;

	mv=matvec(l,u,i,ar,br)-matvec(l,u,i,ai,bi);
	*ri=matvec(l,u,i,ai,br)+matvec(l,u,i,ar,bi);
	*rr=mv;
}
