#include "../real.h"
void newgrn(int n, real_t x[], real_t a[])
{
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int k;

	for (k=n-1; k>=0; k--)
		elmvec(k,n-1,1,a,a,-x[k]);
}
