#include "../real.h"


void lsqortdec(real_t **a, int n, int m, real_t aux[],
					real_t aid[], int ci[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	void ichcol(int, int, int, int, real_t **);
	int j,k,kpiv;
	real_t beta,sigma,norm,w,eps,akk,aidk,temp,*sum;

	sum=allocate_real_vector(1,m);
	norm=0.0;
	aux[3]=m;
	for (k=1; k<=m; k++) {
		w=sum[k]=tammat(1,n,k,k,a,a);
		if (w > norm) norm=w;
	}
	w=aux[5]=sqrt(norm);
	eps=aux[2]*w;
	for (k=1; k<=m; k++) {
		sigma=sum[k];
		kpiv=k;
		for (j=k+1; j<=m; j++)
			if (sum[j] > sigma) {
				sigma=sum[j];
				kpiv=j;
			}
		if (kpiv != k) {
			sum[kpiv]=sum[k];
			ichcol(1,n,k,kpiv,a);
		}
		ci[k]=kpiv;
		akk=a[k][k];
		sigma=tammat(k,n,k,k,a,a);
		w=sqrt(sigma);
		aidk=aid[k]=((akk < 0.0) ? w : -w);
		if (w < eps) {
			aux[3]=k-1;
			break;
		}
		beta=1.0/(sigma-akk*aidk);
		a[k][k]=akk-aidk;
		for (j=k+1; j<=m; j++) {
			elmcol(k,n,j,k,a,a,-beta*tammat(k,n,k,j,a,a));
			temp=a[k][j];
			sum[j] -= temp*temp;
		}
	}
	free_real_vector(sum,1);
}
