#include "../real.h"
int qrisngvaldec(real_t **a, int m, int n, real_t val[], real_t **v,
						real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshreabid(real_t **, int, int, real_t [], real_t [], real_t []);
	void psttfmmat(real_t **, int, real_t **, real_t []);
	void pretfmmat(real_t **, int, int, real_t []);
	int qrisngvaldecbid(real_t [], real_t [], int, int, real_t **,
							real_t **, real_t []);
	int i;
	real_t *b;

	b=allocate_real_vector(1,n);
	hshreabid(a,m,n,val,b,em);
	psttfmmat(a,n,v,b);
	pretfmmat(a,m,n,val);
	i=qrisngvaldecbid(val,b,m,n,a,v,em);
	free_real_vector(b,1);
	return i;
}
