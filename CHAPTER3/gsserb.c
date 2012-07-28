#include "../real.h"
void gsserb(real_t **a, int n, real_t aux[], int ri[], int ci[])
{
	void gsselm(real_t **, int, real_t [], int [], int []);
	real_t onenrminv(real_t **, int);
	void erbelm(int, real_t [], real_t);

	gsselm(a,n,aux,ri,ci);
	if (aux[3] == n) erbelm(n,aux,onenrminv(a,n));
}
