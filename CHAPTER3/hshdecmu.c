#include "../real.h"


void hshdecmul(int n, real_t **a, real_t **b, real_t dwarf)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void hshvecmat(int, int, int, int, real_t, real_t [], real_t **);
	int j,k,k1,n1;
	real_t r,t,c,*v;

	v=allocate_real_vector(1,n);
	k=1;
	n1=n+1;
	for (k1=2; k1<=n1; k1++) {
		r=tammat(k1,n,k,k,b,b);
		if (r > dwarf) {
			r = (b[k][k] < 0.0) ? -sqrt(r+b[k][k]*b[k][k]) :
							sqrt(r+b[k][k]*b[k][k]);
			t=b[k][k]+r;
			c = -t/r;
			b[k][k] = -r;
			v[k]=1.0;
			for (j=k1; j<=n; j++) v[j]=b[j][k]/t;
			hshvecmat(k,n,k1,n,c,v,b);
			hshvecmat(k,n,1,n,c,v,a);
		}
		k=k1;
	}
	free_real_vector(v,1);
}
