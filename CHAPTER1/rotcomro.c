#include "../real.h"
void rotcomrow(int l, int u, int i, int j, real_t **ar, real_t **ai,
					real_t cr, real_t ci, real_t s)
{
	real_t aril,aiil,arjl,aijl;

	for (; l<=u; l++) {
		aril=ar[i][l];
		aiil=ai[i][l];
		arjl=ar[j][l];
		aijl=ai[j][l];
		ar[i][l]=cr*aril+ci*aiil+s*arjl;
		ai[i][l]=cr*aiil-ci*aril+s*aijl;
		ar[j][l]=cr*arjl-ci*aijl-s*aril;
		ai[j][l]=cr*aijl+ci*arjl-s*aiil;
	}
}
