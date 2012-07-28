#include "../real.h"
void chlinv2(real_t **a, int n)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void dupvecrow(int, int, int, real_t [], real_t **);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int i,j,i1;
	real_t r,*u;

	u=allocate_real_vector(1,n);
	for (i=n; i>=1; i--) {
		r=1.0/a[i][i];
		i1=i+1;
		dupvecrow(i1,n,i,u,a);
		for (j=n; j>=i1; j--)
			a[i][j] = -(tamvec(i1,j,j,a,u)+matvec(j+1,n,j,a,u))*r;
		a[i][i]=(r-matvec(i1,n,i,a,u))*r;
	}
	free_real_vector(u,1);
}
