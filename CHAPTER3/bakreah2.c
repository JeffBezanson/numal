#include "../real.h"
void bakreahes2(real_t **a, int n, int n1, int n2, int index[],
						real_t **vec)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t tamvec(int, int, int, real_t **, real_t []);
	void ichrow(int, int, int, int, real_t **);
	int i,l,k;
	real_t *u;

	u=allocate_real_vector(1,n);
	for (i=n; i>=2; i--) {
		for (k=i-2; k>=1; k--) u[k+1]=a[i][k];
		for (k=n1; k<=n2; k++) vec[i][k] += tamvec(2,i-1,k,vec,u);
		l=index[i];
		if (l > i) ichrow(n1,n2,i,l,vec);
	}
	free_real_vector(u,1);
}
