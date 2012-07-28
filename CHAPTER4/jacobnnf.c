#include "../real.h"
void jacobnnf(int n, real_t x[], real_t f[], real_t **jac,
				real_t (*di)(int), void (*funct)(int, real_t[], real_t[]))
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int i,j;
	real_t step,aid,*f1;

	f1=allocate_real_vector(1,n);
	for (i=1; i<=n; i++) {
		step=(*di)(i);
		aid=x[i];
		x[i]=aid+step;
		step=1.0/step;
		(*funct)(n,x,f1);
		for (j=1; j<=n; j++) jac[j][i]=(f1[j]-f[j])*step;
		x[i]=aid;
	}
	free_real_vector(f1,1);
}
