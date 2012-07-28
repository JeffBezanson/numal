#include "../real.h"
void rnk1upd(real_t h[], int n, real_t v[], real_t c)
{
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int j,k;

	k=0;
	j=1;
	do {
		k++;
		elmvec(j,j+k-1,1-j,h,v,v[k]*c);
		j += k;
	} while (k < n);
}
