#include "../real.h"


int hshcomcol(int l, int u, int j, real_t **ar, real_t **ai, real_t tol,
					real_t *k, real_t *c, real_t *s, real_t *t)
{
	void carpol(real_t, real_t, real_t *, real_t *, real_t *);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t vr, mod, h, arlj, ailj;

	vr=tammat(l+1,u,j,j,ar,ar)+tammat(l+1,u,j,j,ai,ai);
	arlj=ar[l][j];
	ailj=ai[l][j];
	carpol(arlj,ailj,&mod,c,s);
	if (vr > tol) {
		vr += arlj*arlj+ailj*ailj;
		h = *k = sqrt(vr);
		*t=vr+mod*h;
		if (arlj == 0.0 && ailj == 0.0)
			ar[l][j]=h;
		else {
			ar[l][j]=arlj + *c * *k;
			ai[l][j]=ailj + *s * *k;
			*s = - *s;
		}
		*c = - *c;
		return (1);
	} else {
		*k=mod;
		*t = -1.0;
		return (0);
	}
}
