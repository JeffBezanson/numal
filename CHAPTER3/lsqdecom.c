#include "../real.h"


void lsqdecomp(real_t **a, int n, int m, int n1, real_t aux[],
					real_t aid[], int ci[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matmat(int, int, int, int, real_t **, real_t **);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	void elmcol(int, int, int, int, real_t **, real_t **, real_t);
	void ichcol(int, int, int, int, real_t **);
	int j,k,kpiv,nr,s,fsum;
	real_t beta,sigma,norm,aidk,akk,w,eps,temp,*sum;

	sum=allocate_real_vector(1,m);
	norm=0.0;
	aux[3]=m;
	nr=n1;
	fsum=1;
	for (k=1; k<=m; k++) {
		if (k == n1+1) {
			fsum=1;
			nr=n;
		}
		if (fsum)
			for (j=k; j<=m; j++) {
				w=sum[j]=tammat(k,nr,j,j,a,a);
				if (w > norm) norm=w;
			}
		fsum=0;
		eps=aux[2]*sqrt(norm);
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
		sigma=tammat(k,nr,k,k,a,a);
		w=sqrt(sigma);
		aidk=aid[k]=((akk < 0.0) ? w : -w);
		if (w < eps) {
			aux[3]=k-1;
			break;
		}
		beta=1.0/(sigma-akk*aidk);
		a[k][k]=akk-aidk;
		for (j=k+1; j<=m; j++) {
			elmcol(k,nr,j,k,a,a,-beta*tammat(k,nr,k,j,a,a));
			temp=a[k][j];
			sum[j] -= temp*temp;
		}
		if (k == n1)
			for (j=n1+1; j<=n; j++)
				for (s=1; s<=m; s++) {
					nr = (s > n1) ? n1 : s-1;
					w=a[j][s]-matmat(1,nr,j,s,a,a);
					a[j][s] = (s > n1) ? w : w/aid[s];
				}
	}
	free_real_vector(sum,1);
}
