#include "../real.h"
int homsol(real_t **a, int m, int n, real_t **v, real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int qrisngvaldec(real_t **, int, int, real_t [], real_t **, real_t []);
	void homsolsvd(real_t **, real_t [], real_t **, int, int);
	int i;
	real_t *val;

	val=allocate_real_vector(1,n);
	i=qrisngvaldec(a,m,n,val,v,em);
	if (i == 0) homsolsvd(a,val,v,m,n);
	free_real_vector(val,1);
	return i;
}

