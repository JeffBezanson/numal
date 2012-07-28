#include "../real.h"
void rotcomcol(int l, int u, int i, int j, real_t **ar, real_t **ai,
					real_t cr, real_t ci, real_t s)
{
	real_t arli,aili,arlj,ailj;

	for (; l<=u; l++) {
		arli=ar[l][i];
		aili=ai[l][i];
		arlj=ar[l][j];
		ailj=ai[l][j];
		ar[l][i]=cr*arli+ci*aili-s*arlj;
		ai[l][i]=cr*aili-ci*arli-s*ailj;
		ar[l][j]=cr*arlj-ci*ailj+s*arli;
		ai[l][j]=cr*ailj+ci*arlj+s*aili;
	}
}
