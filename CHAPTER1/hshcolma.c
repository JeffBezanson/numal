#include "../real.h"
void hshcolmat(int lr, int ur, int lc, int uc, int i,
					real_t x, real_t **u, real_t **a)
{
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	real_t tammat(int, int, int, int, real_t **, real_t **);

	for (; lc<=uc; lc++) elmcol(lr,ur,lc,i,a,u,tammat(lr,ur,lc,i,a,u)*x);
}
