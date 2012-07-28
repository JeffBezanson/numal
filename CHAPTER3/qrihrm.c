#include "../real.h"
int qrihrm(real_t **a, int n, real_t val[], real_t **vr, real_t **vi,
				real_t em[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void hshhrmtri(real_t **, int, real_t [], real_t [], real_t [],
						real_t [], real_t [], real_t []);
	int qrisymtri(real_t **, int, real_t [], real_t [], real_t [], real_t []);
	void bakhrmtri(real_t **, int, int, int, real_t **,
						real_t **, real_t [], real_t []);
	int i,j;
	real_t *b,*bb,*tr,*ti;

	b=allocate_real_vector(1,n);
	bb=allocate_real_vector(1,n);
	tr=allocate_real_vector(1,n-1);
	ti=allocate_real_vector(1,n-1);
	hshhrmtri(a,n,val,b,bb,em,tr,ti);
	for (i=1; i<=n; i++) {
		vr[i][i]=1.0;
		for (j=i+1; j<=n; j++) vr[i][j]=vr[j][i]=0.0;
	}
	b[n]=bb[n]=0.0;
	i=qrisymtri(vr,n,val,b,bb,em);
	bakhrmtri(a,n,i+1,n,vr,vi,tr,ti);
	free_real_vector(b,1);
	free_real_vector(bb,1);
	free_real_vector(tr,1);
	free_real_vector(ti,1);
	return i;
}

