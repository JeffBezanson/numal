#include "../real.h"
int psdinv(real_t **a, int m, int n, real_t em[])
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	int qrisngvaldec(real_t **, int, int, real_t [], real_t **, real_t []);
	void psdinvsvd(real_t **, real_t [], real_t **, int, int, real_t []);
	int i;
	real_t *val,**v;

	val=allocate_real_vector(1,n);
	v=allocate_real_matrix(1,n,1,n);
	i=qrisngvaldec(a,m,n,val,v,em);
	if (i == 0) psdinvsvd(a,val,v,m,n,em);
	free_real_vector(val,1);
	free_real_matrix(v,1,n,1);
	return i;
}

