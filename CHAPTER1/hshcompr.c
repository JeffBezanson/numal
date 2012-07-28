#include "../real.h"
void hshcomprd(int i, int ii, int l, int u, int j, real_t **ar,
					real_t **ai,	real_t **br, real_t **bi, real_t t)
{
	void elmcomcol(int, int, int, int, real_t **, real_t **,
						real_t **, real_t **, real_t, real_t);
	real_t tammat(int, int, int, int, real_t **, real_t **);

	for (; l<=u; l++)
		elmcomcol(i,ii,l,j,ar,ai,br,bi,
					(-tammat(i,ii,j,l,br,ar)-tammat(i,ii,j,l,bi,ai))/t,
					(tammat(i,ii,j,l,bi,ar)-tammat(i,ii,j,l,br,ai))/t);
}
