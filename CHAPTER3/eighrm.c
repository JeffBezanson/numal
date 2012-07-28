#include "../real.h"
void eighrm(real_t **a, int n, int numval, real_t val[], real_t **vecr,
				real_t **veci, real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshhrmtri(real_t **, int, real_t [], real_t [], real_t [],
						real_t [], real_t [], real_t []);
	void valsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t []);
	void vecsymtri(real_t [], real_t [], int, int, int,
						real_t [], real_t **, real_t []);
	void bakhrmtri(real_t **, int, int, int, real_t **,
						real_t **, real_t [], real_t []);
	real_t *bb,*tr,*ti,*d,*b;

	bb=allocate_real_vector(1,n-1);
	tr=allocate_real_vector(1,n-1);
	ti=allocate_real_vector(1,n-1);
	d=allocate_real_vector(1,n);
	b=allocate_real_vector(1,n);
	hshhrmtri(a,n,d,b,bb,em,tr,ti);
	valsymtri(d,bb,n,1,numval,val,em);
	b[n]=0.0;
	vecsymtri(d,b,n,1,numval,val,vecr,em);
	bakhrmtri(a,n,1,numval,vecr,veci,tr,ti);
	free_real_vector(bb,1);
	free_real_vector(tr,1);
	free_real_vector(ti,1);
	free_real_vector(d,1);
	free_real_vector(b,1);
}

