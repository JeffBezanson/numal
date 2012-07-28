#include "../real.h"
void homsolsvd(real_t **u, real_t val[], real_t **v, int m, int n)
{
	void ichcol(int, int, int, int, real_t **);
	int i,j;
	real_t x;

	for (i=n; i>=2; i--)
		for (j=i-1; j>=1; j--)
			if (val[i] > val[j]) {
				x=val[i];
				val[i]=val[j];
				val[j]=x;
				ichcol(1,m,i,j,u);
				ichcol(1,n,i,j,v);
			}
}

