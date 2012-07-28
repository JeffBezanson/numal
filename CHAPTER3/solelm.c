#include "../real.h"
void solelm(real_t **a, int n, int ri[], int ci[], real_t b[])
{
	void sol(real_t **, int, int [], real_t []);
	int r,cir;
	real_t w;

	sol(a,n,ri,b);
	for (r=n; r>=1; r--) {
		cir=ci[r];
		if (cir != r) {
			w=b[r];
			b[r]=b[cir];
			b[cir]=w;
		}
	}
}
