#include "../real.h"
void hshvectam(int lr, int ur, int lc, int uc,
					real_t x, real_t u[], real_t **a)
{
	real_t matvec(int, int, int, real_t **, real_t []);
	void elmrowvec(int, int, int, real_t **, real_t [], real_t);

	for (; lr<=ur; lr++) elmrowvec(lc,uc,lr,a,u,matvec(lc,uc,lr,a,u)*x);
}
