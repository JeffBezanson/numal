#include "../real.h"
void gssitisolerb(real_t **a, int n, real_t aux[], real_t b[])
{
	int *allocate_integer_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_matrix(real_t **, int, int, int);
	void gssnri(real_t **, int, real_t [], int [], int []);
	void itisolerb(real_t **, real_t **, int, real_t [],
						int [], int [], real_t []);
	void dupmat(int, int, int, int, real_t **, real_t **);
	int *ri,*ci;
	real_t **aa;

	ri=allocate_integer_vector(1,n);
	ci=allocate_integer_vector(1,n);
	aa=allocate_real_matrix(1,n,1,n);

	dupmat(1,n,1,n,aa,a);
	gssnri(a,n,aux,ri,ci);
	if (aux[3] == n) itisolerb(aa,a,n,aux,ri,ci,b);

	free_integer_vector(ri,1);
	free_integer_vector(ci,1);
	free_real_matrix(aa,1,n,1);
}
