#include "../real.h"
void lsqortdecsol(real_t **a, int n, int m, real_t aux[],
						real_t diag[], real_t b[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void lsqortdec(real_t **, int, int, real_t [], real_t [], int []);
	void lsqdglinv(real_t **, int, real_t [], int [], real_t []);
	void lsqsol(real_t **, int, int, real_t [], int [], real_t []);
	int *ci;
	real_t *aid;

	ci=allocate_integer_vector(1,m);
	aid=allocate_real_vector(1,m);
	lsqortdec(a,n,m,aux,aid,ci);
	if (aux[3] == m) {
		lsqdglinv(a,m,aid,ci,diag);
		lsqsol(a,n,m,aid,ci,b);
	}
	free_integer_vector(ci,1);
	free_real_vector(aid,1);
}
