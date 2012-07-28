#include "../real.h"
void gssnri(real_t **a, int n, real_t aux[], int ri[], int ci[])
{
	void gsselm(real_t **, int, real_t [], int [], int []);
	real_t onenrminv(real_t **, int);

	gsselm(a,n,aux,ri,ci);
	if (aux[3] == n) aux[9]=onenrminv(a,n);
}
