#include "../real.h"


int reaeigval(real_t **a, int n, real_t em[], real_t val[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void tfmreahes(real_t **, int, real_t [], int []);
	void eqilbr(real_t **, int, real_t [], real_t [], int []);
	int reavalqri(real_t **, int, real_t [], real_t []);
	int i,j,k,*ind,*ind0;
	real_t r,*d;

	ind=allocate_integer_vector(1,n);
	ind0=allocate_integer_vector(1,n);
	d=allocate_real_vector(1,n);
	eqilbr(a,n,em,d,ind0);
	tfmreahes(a,n,em,ind);
	k=reavalqri(a,n,em,val);
	for (i=k+1; i<=n; i++)
		for (j=i+1; j<=n; j++)
			if (val[j] > val[i]) {
				r=val[i];
				val[i]=val[j];
				val[j]=r;
			}
	free_integer_vector(ind,1);
	free_integer_vector(ind0,1);
	free_real_vector(d,1);
	return k;
}

