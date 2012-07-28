#include "../real.h"
void solsvdovr(real_t **u, real_t val[], real_t **v, int m, int n,
					real_t x[], real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int i;
	real_t min,*x1;

	x1=allocate_real_vector(1,n);
	min=em[6];
	for (i=1; i<=n; i++)
		x1[i] = (val[i] <= min) ? 0.0 : tamvec(1,m,i,u,x)/val[i];
	for (i=1; i<=n; i++) x[i]=matvec(1,n,i,v,x1);
	free_real_vector(x1,1);
}
