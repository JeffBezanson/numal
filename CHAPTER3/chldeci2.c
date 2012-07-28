#include "../real.h"
void chldecinv2(real_t **a, int n, real_t aux[])
{
	void chldec2(real_t **, int, real_t []);
	void chlinv2(real_t **, int);

	chldec2(a,n,aux);
	if (aux[3] == n) chlinv2(a,n);
}
