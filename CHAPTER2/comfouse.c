#include "../real.h"


void comfouser(int n, real_t theta, real_t a[], real_t *rr, real_t *ri)
{
	int k;
	real_t c,cc,lambda,h,dun,un,un1,temp;

	c=cos(theta);
	if (c < -0.5) {
		temp=cos(theta/2.0);
		lambda=4.0*temp*temp;
		un=dun=0.0;
		for (k=n; k>=0; k--) {
			un=dun-un;
			dun=lambda*un-dun+a[k];
		}
		*rr=dun-lambda/2.0*un;
	} else {
		if (c > 0.5) {
			temp=sin(theta/2.0);
			lambda = -4.0*temp*temp;
			un=dun=0.0;
			for (k=n; k>=0; k--) {
				un += dun;
				dun += lambda*un+a[k];
			}
			*rr=dun-lambda/2.0*un;
		} else {
			cc=c+c;
			un=un1=0.0;
			for (k=n; k>=1; k--) {
				h=cc*un-un1+a[k];
				un1=un;
				un=h;
			}
			*rr=a[0]+un*c-un1;
		}
	}
	*ri=un*sin(theta);
}
