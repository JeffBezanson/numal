#include "../real.h"


void symeigimp(int n, real_t **a, real_t **vec, real_t val[],
					real_t lbound[], real_t ubound[], real_t aux[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	void orthog(int, int, int, real_t **);
	void mergesort(real_t [], int [], int, int);
	void vecperm(int [], int, int, real_t []);
	void rowperm(int [], int, int, int, real_t **);
	int qrisym(real_t **, int, real_t [], real_t []);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t infnrmmat(int, int, int, int, int *, real_t **);
	int k,i,j,i0,i1,i01,iter,maxitp1,stop,n1,i0m1,i1p1,*perm;
	real_t s,max,tol,mateps,relerra,reltolr,norma,eps2,dl,dr,m1,
			**r,**p,**y,**pp,*rq,*eps,*z,*val3,*val4,*eta,em[6];
	double dtemp;

	perm=allocate_integer_vector(1,n);
	rq=allocate_real_vector(1,n);
	eps=allocate_real_vector(1,n);
	z=allocate_real_vector(1,n);
	val3=allocate_real_vector(1,n);
	eta=allocate_real_vector(1,n);
	r=allocate_real_matrix(1,n,1,n);
	p=allocate_real_matrix(1,n,1,n);
	y=allocate_real_matrix(1,n,1,n);

	norma=infnrmmat(1,n,1,n,&i,a);
	relerra=aux[0];
	reltolr=aux[2];
	maxitp1=aux[4]+1.0;
	mateps=relerra*norma;
	tol=reltolr*norma;
	for (iter=1; iter<=maxitp1; iter++) {
		if (iter == 1)
			stop=0;
		else
			stop=1;
		max=0.0;
		for (j=1; j<=n; j++)
			for (i=1; i<=n; i++) {
				dtemp = -(double)(vec[i][j])*(double)(val[j]);
				for (k=1; k<=n; k++)
					dtemp += (double)(a[i][k])*(double)(vec[k][j]);
				r[i][j]=dtemp;
				if (fabs(r[i][j]) > max) max=fabs(r[i][j]);
			}
		if (max > tol) stop=0;
		if ((!stop) && (iter < maxitp1)) {
			for (i=1; i<=n; i++) {
				dtemp=(double)(val[i]);
				for (k=1; k<=n; k++)
					dtemp += (double)(vec[k][i])*(double)(r[k][i]);
				rq[i]=dtemp;
			}
			for (j=1; j<=n; j++) {
				for (i=1; i<=n; i++)
					eta[i]=r[i][j]-(rq[j]-val[j])*vec[i][j];
				z[j]=sqrt(vecvec(1,n,0,eta,eta));
			}
			mergesort(rq,perm,1,n);
			vecperm(perm,1,n,rq);
			for (i=1; i<=n; i++) {
				eps[i]=z[perm[i]];
				val3[i]=val[perm[i]];
				rowperm(perm,1,n,i,vec);
				rowperm(perm,1,n,i,r);
			}
			for (i=1; i<=n; i++)
				for (j=i; j<=n; j++)
					p[i][j]=p[j][i]=tammat(1,n,i,j,vec,r);
		}
		i0=1;
		do {
			j=i1=i0;
			j++;
			while ((j > n) ? 0 :
						(rq[j]-rq[j-1] <= sqrt((eps[j]+eps[j-1])*norma))) {
				i1=j;
				j++;
			}
			if (stop || (iter == maxitp1)) {
				i=i0;
				do {
					j=i01=i;
					j++;
					while ((j>i1) ? 0 : rq[j]-rq[j-1] <= eps[j]+eps[j-1]) {
						i01=j;
						j++;
					}
					if (i == i01) {
						if (i < n) {
							if (i == 1)
								dl=dr=rq[i+1]-rq[i]-eps[i+1];
							else {
								dl=rq[i]-rq[i-1]-eps[i-1];
								dr=rq[i+1]-rq[i]-eps[i+1];
							}
						} else
							dl=dr=rq[i]-rq[i-1]-eps[i-1];
						eps2=eps[i]*eps[i];
						lbound[i]=eps2/dr+mateps;
						ubound[i]=eps2/dl+mateps;
					} else
						for (k=i; k<=i01; k++)
							lbound[k]=ubound[k]=eps[k]+mateps;
					i01++;
					i=i01;
				} while (i <= i1);	/* bounds */
			} else {
				if (i0 == i1) {
					for (k=1; k<=n; k++)
						if (k == i0)
							y[k][i0]=1.0;
						else
							r[k][i0]=p[k][i0];
					val[i0]=rq[i0];
				} else {
					n1=i1-i0+1;
					em[0]=em[2]=FLT_EPSILON;
					em[4]=10*n1;
					val4=allocate_real_vector(1,n1);
					pp=allocate_real_matrix(1,n1,1,n1);
					m1=0.0;
					for (k=i0; k<=i1; k++) m1 += val3[k];
					m1 /= n1;
					for (i=1; i<=n1; i++)
						for (j=1; j<=n1; j++) {
							pp[i][j]=p[i+i0-1][j+i0-1];
							if (i == j) pp[i][j] += val3[j+i0-1]-m1;
						}
					for (i=i0; i<=i1; i++) {
						val3[i]=m1;
						val[i]=rq[i];
					}
					qrisym(pp,n1,val4,em);
					mergesort(val4,perm,1,n1);
					for (i=1; i<=n1; i++)
						for (j=1; j<=n1; j++)
							p[i+i0-1][j+i0-1]=pp[i][perm[j]];
					i0m1=i0-1;
					i1p1=i1+1;
					for (j=i0; j<=i1; j++) {
						for (i=1; i<=i0m1; i++) {
							s=0.0;
							for (k=i0; k<=i1; k++) s += p[i][k]*p[k][j];
							r[i][j]=s;
						}
						for (i=i1p1; i<=n; i++) {
							s=0.0;
							for (k=i0; k<=i1; k++) s += p[i][k]*p[k][j];
							r[i][j]=s;
						}
						for (i=i0; i<=i1; i++) y[i][j]=p[i][j];
					}
					free_real_vector(val4,1);
					free_real_matrix(pp,1,n1,1);
				}	/* innerblock */
			}	/* not stop */
			i0=i1+1;
		} while (i0 <= n);	/* while i0 loop */
		if ((!stop) && (iter < maxitp1)) {
			for (j=1; j<=n; j++)
				for (i=1; i<=n; i++)
					if (val3[i] != val3[j])
						y[i][j]=r[i][j]/(val3[j]-val3[i]);
			for (i=1; i<=n; i++) {
				for (j=1; j<=n; j++) z[j]=matmat(1,n,i,j,vec,y);
				for (j=1; j<=n; j++) vec[i][j]=z[j];
			}
			orthog(n,1,n,vec);
		} else {
			aux[5]=iter-1;
			break;
		}
	}	/* for iter loop */
	aux[1]=norma;
	aux[3]=max;

	free_integer_vector(perm,1);
	free_real_vector(rq,1);
	free_real_vector(eps,1);
	free_real_vector(z,1);
	free_real_vector(val3,1);
	free_real_vector(eta,1);
	free_real_matrix(r,1,n,1);
	free_real_matrix(p,1,n,1);
	free_real_matrix(y,1,n,1);
}

