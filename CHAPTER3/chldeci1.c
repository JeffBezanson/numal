#include "../real.h"
void chldecinv1(real_t a[], int n, real_t aux[])
{
	void chldec1(real_t [], int, real_t []);
	void chlinv1(real_t [], int);

	chldec1(a,n,aux);
	if (aux[3] == n) chlinv1(a,n);
}
