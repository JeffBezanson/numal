#include "../real.h"


real_t sinser(int n, real_t theta, real_t b[])
{
	int k;
	real_t c,cc,lambda,h,dun,un,un1,temp;

	c=cos(theta);
	if (c < -0.5) {
		temp=cos(theta/2.0);
		lambda=4.0*temp*temp;
		un=dun=0.0;
		for (k=n; k>=1; k--) {
			dun=lambda*un-dun+b[k];
			un=dun-un;
		}
	} else {
		if (c > 0.5) {
			temp=sin(theta/2.0);
			lambda = -4.0*temp*temp;
			un=dun=0.0;
			for (k=n; k>=1; k--) {
				dun += lambda*un+b[k];
				un += dun;
			}
		} else {
			cc=c+c;
			un=un1=0.0;
			for (k=n; k>=1; k--) {
				h=cc*un-un1+b[k];
				un1=un;
				un=h;
			}
		}
	}
	return (un*sin(theta));
}
