#include "../real.h"
void soltripiv(real_t sub[], real_t diag[], real_t super[], int n,
					real_t aid[], int piv[], real_t b[])
{
	int i,n1;
	real_t bi,bi1,r,s,t;

	n1=n-1;
	for (i=1; i<=n1; i++) {
		if (piv[i]) {
			bi=b[i+1];
			bi1=b[i];
		} else {
			bi=b[i];
			bi1=b[i+1];
		}
		r=b[i]=bi/diag[i];
		b[i+1]=bi1-sub[i]*r;
	}
	r = b[n] /= diag[n];
	t = b[n1] -= super[n1]*r;
	for (i=n-2; i>=1; i--) {
		s=r;
		r=t;
		t = b[i] -= super[i]*r + ((piv[i]) ? aid[i]*s : 0.0);
	}
}
