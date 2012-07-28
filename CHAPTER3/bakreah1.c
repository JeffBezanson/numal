#include "../real.h"
void bakreahes1(real_t **a, int n, int index[], real_t v[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matvec(int, int, int, real_t **, real_t []);
	int i,l;
	real_t w,*x;

	x=allocate_real_vector(1,n);
	for (i=2; i<=n; i++) x[i-1]=v[i];
	for (i=n; i>=2; i--) {
		v[i] += matvec(1,i-2,i,a,x);
		l=index[i];
		if (l > i) {
			w=v[i];
			v[i]=v[l];
			v[l]=w;
		}
	}
	free_real_vector(x,1);
}
