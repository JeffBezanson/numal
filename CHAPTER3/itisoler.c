#include "../real.h"


void itisolerb(real_t **a, real_t **lu, int n, real_t aux[],
				int ri[], int ci[], real_t b[])
{
	void itisol(real_t **, real_t **, int, real_t [],
					int [], int [], real_t []);
	int i;
	real_t nrmsol,nrminv,nrmb,alfa,tola,eps;

	eps=aux[0];
	nrminv=aux[9];
	tola=aux[5]*aux[6];
	nrmb=nrmsol=0.0;
	for (i=1; i<=n; i++) nrmb += fabs(b[i]);
	itisol(a,lu,n,aux,ri,ci,b);
	for (i=1; i<=n; i++) nrmsol += fabs(b[i]);
	alfa=1.0-(1.06*eps*aux[7]*(0.75*n+4.5)*n*n+tola)*nrminv;
	if (alfa < eps)
		aux[11] = -1.0;
	else {
		alfa=((aux[13]+aux[8]*nrmb)/nrmsol+tola)*nrminv/alfa;
		aux[11]=(1.0-alfa < eps) ? -1.0 : alfa/(1.0-alfa);
	}
}
