#include "../real.h"
void baklbr(int n, int n1, int n2, real_t d[], int inter[], real_t **vec)
{
	void ichrow(int, int, int, int, real_t **);
	int i,j,k,p,q;
	real_t di;

	p=1;
	q=n;
	for (i=1; i<=n; i++) {
		di=d[i];
		if (di != 1)
			for (j=n1; j<=n2; j++) vec[i][j] *= di;
		k=inter[i];
		if (k > 0)
			p++;
		else
			if (k < 0) q--;
	}
	for (i=p-1+n-q; i>=1; i--) {
		k=inter[i];
		if (k > 0) {
			p--;
			if (k != p) ichrow(n1,n2,k,p,vec);
		} else {
			q++;
			if (-k != q) ichrow(n1,n2,-k,q,vec);
		}
	}
}
