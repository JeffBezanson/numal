#include "../real.h"
void gsssol(real_t **a, int n, real_t aux[], real_t b[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void solelm(real_t **, int, int [], int [], real_t []);
	void gsselm(real_t **, int, real_t [], int [], int []);
	int *ri,*ci;

	ri=allocate_integer_vector(1,n);
	ci=allocate_integer_vector(1,n);
	gsselm(a,n,aux,ri,ci);
	if (aux[3] == n) solelm(a,n,ri,ci,b);
	free_integer_vector(ri,1);
	free_integer_vector(ci,1);
}
