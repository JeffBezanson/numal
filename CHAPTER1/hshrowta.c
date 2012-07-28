#include "../real.h"
void hshrowtam(int lr, int ur, int lc, int uc, int i,
					real_t x, real_t **u, real_t **a)
{
	real_t mattam(int, int, int, int, real_t **, real_t **);
	void elmrow(int, int, int, int, real_t **, real_t **, real_t);

	for (; lr<=ur; lr++) elmrow(lc,uc,lr,i,a,u,mattam(lc,uc,lr,i,a,u)*x);
}
