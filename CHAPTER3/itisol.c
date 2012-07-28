#include "../real.h"


void itisol(real_t **a, real_t **lu, int n, real_t aux[],
				int ri[], int ci[], real_t b[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void solelm(real_t **, int, int [], int [], real_t []);
	void inivec(int, int, real_t [], real_t);
	void dupvec(int, int, int, real_t [], real_t []);
	int i,j,iter,maxiter;
	real_t maxerx,erx,nrmres,nrmsol,r,rr,*res,*sol;
	double dtemp;

	res=allocate_real_vector(1,n);
	sol=allocate_real_vector(1,n);
	maxerx=erx=aux[10];
	maxiter=aux[12];
	inivec(1,n,sol,0.0);
	dupvec(1,n,0,res,b);
	iter=1;
	do {
		solelm(lu,n,ri,ci,res);
		erx=nrmsol=nrmres=0.0;
		for (i=1; i<=n; i++) {
			r=res[i];
			erx += fabs(r);
			rr=sol[i]+r;
			sol[i]=rr;
			nrmsol += fabs(rr);
		}
		erx /= nrmsol;
		for (i=1; i<=n; i++) {
			dtemp = -(double)b[i];
			for (j=1; j<=n; j++)
				dtemp += (double)a[i][j]*(double)sol[j];
			r = -dtemp;
			res[i]=r;
			nrmres += fabs(r);
		}
		iter++;
	} while ((iter <= maxiter) && (maxerx < erx));
	dupvec(1,n,0,b,sol);
	aux[11]=erx;
	aux[13]=nrmres;
	free_real_vector(res,1);
	free_real_vector(sol,1);
}
