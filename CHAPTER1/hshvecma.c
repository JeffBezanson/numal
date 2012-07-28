#include "../real.h"
void hshvecmat(int lr, int ur, int lc, int uc,
					real_t x, real_t u[], real_t **a)
{
	void elmcolvec(int, int, int, real_t **, real_t [], real_t);
	real_t tamvec(int, int, int, real_t **, real_t []);

	for (; lc<=uc; lc++) elmcolvec(lr,ur,lc,a,u,tamvec(lr,ur,lc,a,u)*x);
}
