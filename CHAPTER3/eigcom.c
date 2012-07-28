#include "../real.h"
int eigcom(real_t **ar, real_t **ai, int n, real_t em[], real_t valr[],
				real_t vali[], real_t **vr, real_t **vi)
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void eqilbrcom(real_t **, real_t **, int, real_t [], real_t [], int []);
	real_t comeucnrm(real_t **, real_t **, int, int);
	void hshcomhes(real_t **, real_t **, int, real_t [], real_t [],
						real_t [], real_t [], real_t []);
	int qricom(real_t **, real_t **, real_t [], int, real_t [],
					real_t [], real_t [], real_t **, real_t **);
	void bakcomhes(real_t **, real_t **, real_t [], real_t [], real_t [],
						real_t **, real_t **, int, int, int);
	void baklbrcom(int, int, int, real_t [], int [], real_t **, real_t **);
	void sclcom(real_t **, real_t **, int, int, int);
	int i,*ind;
	real_t *d,*b,*del,*tr,*ti;

	ind=allocate_integer_vector(1,n);
	d=allocate_real_vector(1,n);
	b=allocate_real_vector(1,n);
	del=allocate_real_vector(1,n);
	tr=allocate_real_vector(1,n);
	ti=allocate_real_vector(1,n);
	eqilbrcom(ar,ai,n,em,d,ind);
	em[1]=comeucnrm(ar,ai,n-1,n);
	hshcomhes(ar,ai,n,em,b,tr,ti,del);
	i=qricom(ar,ai,b,n,em,valr,vali,vr,vi);
	if (i == 0) {
		bakcomhes(ar,ai,tr,ti,del,vr,vi,n,1,n);
		baklbrcom(n,1,n,d,ind,vr,vi);
		sclcom(vr,vi,n,1,n);
	}
	free_integer_vector(ind,1);
	free_real_vector(d,1);
	free_real_vector(b,1);
	free_real_vector(del,1);
	free_real_vector(tr,1);
	free_real_vector(ti,1);
	return i;
}

