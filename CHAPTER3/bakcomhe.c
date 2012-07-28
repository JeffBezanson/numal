#include "../real.h"
void bakcomhes(real_t **ar, real_t **ai, real_t tr[], real_t ti[],
					real_t del[], real_t **vr, real_t **vi, int n,
					int n1, int n2)
{
	void hshcomprd(int, int, int, int, int, real_t **,
						real_t **, real_t **, real_t **, real_t);
	void comrowcst(int, int, int, real_t **, real_t **, real_t, real_t);
	int i,r,rm1;
	real_t h;

	for (i=2; i<=n; i++) comrowcst(n1,n2,i,vr,vi,tr[i],ti[i]);
	r=n-1;
	for (rm1=n-2; rm1>=1; rm1--) {
		h=del[rm1];
		if (h > 0.0) hshcomprd(r,n,n1,n2,rm1,vr,vi,ar,ai,h);
		r=rm1;
	}
}
