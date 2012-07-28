#include "../real.h"


void erbelm(int n, real_t aux[], real_t nrminv)
{
	real_t aid,eps;

	eps=aux[0];
	aid=(1.06*eps*(0.75*n+4.5)*(n*n)*aux[7]+aux[5]*aux[6])*nrminv;
	aux[11]=(2.0*aid >= (1.0-eps)) ? -1.0 : aid/(1.0-2.0*aid);
	aux[9]=nrminv;
}
