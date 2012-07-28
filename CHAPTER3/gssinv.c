#include "../real.h"
void gssinv(real_t **a, int n, real_t aux[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void gsselm(real_t **, int, real_t [], int [], int []);
	real_t inv1(real_t **, int, int [], int [], int);
	int *ri,*ci;

	ri=allocate_integer_vector(1,n);
	ci=allocate_integer_vector(1,n);
	gsselm(a,n,aux,ri,ci);
	if (aux[3] == n) aux[9]=inv1(a,n,ri,ci,1);
	free_integer_vector(ri,1);
	free_integer_vector(ci,1);
}
