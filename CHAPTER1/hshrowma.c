#include "../real.h"
void hshrowmat(int lr, int ur, int lc, int uc, int i,
					real_t x, real_t **u, real_t **a)
{
	real_t matmat(int, int, int, int, real_t **, real_t **);
	void elmcolrow(int, int, int, int, real_t **, real_t **, real_t);

	for (; lc<=uc; lc++) elmcolrow(lr,ur,lc,i,a,u,matmat(lr,ur,i,lc,u,a)*x);
}
