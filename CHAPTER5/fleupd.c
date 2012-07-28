#include "../real.h"
void fleupd(real_t h[], int n, real_t v[], real_t w[], real_t c1, real_t c2)
{
	int i,j,k;
	real_t vk,wk;

	k=0;
	j=1;
	do {
		k++;
		vk = -w[k]*c1+v[k]*c2;
		wk=v[k]*c1;
		for (i=0; i<=k-1; i++) h[i+j] += v[i+1]*vk-w[i+1]*wk;
		j += k;
	} while (k < n);
}
