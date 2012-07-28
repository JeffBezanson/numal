#include "../real.h"


void comfouser1(int n, real_t theta, real_t ar[], real_t ai[],
					real_t *rr, real_t *ri)
{
	int k;
	real_t h,hr,hi,co,si;

	hr=hi=0.0;
	co=cos(theta);
	si=sin(theta);
	for (k=n; k>=1; k--) {
		h=co*hr-si*hi+ar[k];
		hi=co*hi+si*hr+ai[k];
		hr=h;
	}
	*rr=co*hr-si*hi+ar[0];
	*ri=co*hi+si*hr+ai[0];
}
