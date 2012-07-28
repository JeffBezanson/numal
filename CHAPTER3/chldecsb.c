#include "../real.h"
void chldecsolbnd(real_t a[], int n, int w, real_t aux[], real_t b[])
{
	void chldecbnd(real_t [], int, int, real_t []);
	void chlsolbnd(real_t [], int, int, real_t []);

	chldecbnd(a,n,w,aux);
	if (aux[3] == n) chlsolbnd(a,n,w,b);
}
