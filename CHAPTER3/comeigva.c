#include "../real.h"
int comeigval(real_t **a, int n, real_t em[], real_t re[], real_t im[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void eqilbr(real_t **, int, real_t [], real_t [], int []);
	void tfmreahes(real_t **, int, real_t [], int []);
	int comvalqri(real_t **, int, real_t [], real_t [], real_t []);
	int i,*ind,*ind0;
	real_t *d;

	ind=allocate_integer_vector(1,n);
	ind0=allocate_integer_vector(1,n);
	d=allocate_real_vector(1,n);
	eqilbr(a,n,em,d,ind0);
	tfmreahes(a,n,em,ind);
	i=comvalqri(a,n,em,re,im);
	free_integer_vector(ind,1);
	free_integer_vector(ind0,1);
	free_real_vector(d,1);
	return i;
}
