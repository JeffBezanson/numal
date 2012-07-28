#include "../real.h"
int reaeig3(real_t **a, int n, real_t em[], real_t val[], real_t **vec)
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void tfmreahes(real_t **, int, real_t [], int []);
	void bakreahes2(real_t **, int, int, int, int [], real_t **);
	void eqilbr(real_t **, int, real_t [], real_t [], int []);
	void baklbr(int, int, int, real_t [], int [], real_t **);
	void reascl(real_t **, int, int, int);
	int reaqri(real_t **, int, real_t [], real_t [], real_t **);
	int i,*ind,*ind0;
	real_t *d;

	ind=allocate_integer_vector(1,n);
	ind0=allocate_integer_vector(1,n);
	d=allocate_real_vector(1,n);
	eqilbr(a,n,em,d,ind0);
	tfmreahes(a,n,em,ind);
	i=reaqri(a,n,em,val,vec);
	if (i == 0) {
		bakreahes2(a,n,1,n,ind,vec);
		baklbr(n,1,n,d,ind0,vec);
		reascl(vec,n,1,n);
	}
	free_integer_vector(ind,1);
	free_integer_vector(ind0,1);
	free_real_vector(d,1);
	return i;
}
