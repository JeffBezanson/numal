#include "../real.h"
void hshcoltam(int lr, int ur, int lc, int uc, int i,
					real_t x, real_t **u, real_t **a)
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void elmrowcol(int, int, int, int, real_t **, real_t **, real_t);

	for (; lr<=ur; lr++) elmrowcol(lc,uc,lr,i,a,u,matmat(lc,uc,lr,i,a,u)*x);
}
