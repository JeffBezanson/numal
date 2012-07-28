#include "../real.h"
void psdinvsvd(real_t **u, real_t val[], real_t **v, int m, int n,
					real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matvec(int, int, int, real_t **, real_t []);
	int i,j;
	real_t min,vali,*x;

	x=allocate_real_vector(1,n);
	min=em[6];
	for (i=1; i<=n; i++)
		if (val[i] > min) {
			vali=1.0/val[i];
			for (j=1; j<=m; j++) u[j][i] *= vali;
		} else
			for (j=1; j<=m; j++) u[j][i]=0.0;
	for (i=1; i<=m; i++) {
		for (j=1; j<=n; j++) x[j]=u[i][j];
		for (j=1; j<=n; j++) u[i][j]=matvec(1,n,j,v,x);
	}
	free_real_vector(x,1);
}

