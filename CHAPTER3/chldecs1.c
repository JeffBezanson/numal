#include "../real.h"
void chldecsol1(real_t a[], int n, real_t aux[], real_t b[])
{
	void chldec1(real_t [], int, real_t []);
	void chlsol1(real_t [], int, real_t []);

	chldec1(a,n,aux);
	if (aux[3] == n) chlsol1(a,n,b);
}
