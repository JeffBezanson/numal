#include "../real.h"
void chlinv1(real_t a[], int n)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t seqvec(int, int, int, int, real_t [], real_t []);
	real_t symmatvec(int, int, int, real_t [], real_t []);
	int i,ii,i1,j,ij;
	real_t r,*u;

	u=allocate_real_vector(1,n);
	ii=((n+1)*n)/2;
	for (i=n; i>=1; i--) {
		r=1.0/a[ii];
		i1=i+1;
		ij=ii+i;
		for (j=i1; j<=n; j++) {
			u[j]=a[ij];
			ij += j;
		}
		for (j=n; j>=i1; j--) {
			ij -= j;
			a[ij] = -symmatvec(i1,n,j,a,u)*r;
		}
		a[ii]=(r-seqvec(i1,n,ii+i,0,a,u))*r;
		ii -= i;
	}
	free_real_vector(u,1);
}
