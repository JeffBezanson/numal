#include "../real.h"
void eigsym2(real_t **a, int n, int numval, real_t val[],
				real_t **vec, real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void tfmsymtri2(real_t **, int, real_t [], real_t [], real_t [],
						real_t []);
	void valsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t []);
	void vecsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t **, real_t []);
	void baksymtri2(real_t **, int, int, int, real_t **);
	real_t *b,*bb,*d;

	b=allocate_real_vector(1,n);
	bb=allocate_real_vector(1,n);
	d=allocate_real_vector(1,n);
	tfmsymtri2(a,n,d,b,bb,em);
	valsymtri(d,bb,n,1,numval,val,em);
	vecsymtri(d,b,n,1,numval,val,vec,em);
	baksymtri2(a,n,1,numval,vec);
	free_real_vector(b,1);
	free_real_vector(bb,1);
	free_real_vector(d,1);
}

