#include "../real.h"
void chldecsol2(real_t **a, int n, real_t aux[], real_t b[])
{
	void chldec2(real_t **, int, real_t []);
	void chlsol2(real_t **, int, real_t []);

	chldec2(a,n,aux);
	if (aux[3] == n) chlsol2(a,n,b);
}
