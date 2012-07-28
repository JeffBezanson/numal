#include "../real.h"
void davupd(real_t h[], int n, real_t v[], real_t w[], real_t c1, real_t c2)
{
	int i,j,k;
	real_t vk,wk;

	k=0;
	j=1;
	do {
		k++;
		vk=v[k]*c1;
		wk=w[k]*c2;
		for (i=0; i<=k-1; i++) h[i+j] += v[i+1]*vk-w[i+1]*wk;
		j += k;
	} while (k < n);
}
