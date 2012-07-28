#include "../real.h"


real_t fouser(int n, real_t theta, real_t a[])
{
	int k;
	real_t c,cc,lambda,h,dun,un,un1,c2,s2;

	c=cos(theta);
	if (c < -0.5) {
		c2=cos(theta/2.0);
		lambda=4.0*c2*c2;
		un=dun=0.0;
		for (k=n; k>=0; k--) {
			un=dun-un;
			dun=lambda*un-dun+a[k];
		}
		return (dun+2.0*c2*(sin(theta/2.0)-c2)*un);
	} else {
		if (c > 0.5) {
			s2=sin(theta/2.0);
			lambda = -4.0*s2*s2;
			un=dun=0.0;
			for (k=n; k>=0; k--) {
				un += dun;
				dun += lambda*un+a[k];
			}
			return (dun+2.0*s2*(s2+cos(theta/2.0))*un);
		} else {
			cc=c+c;
			un=un1=0.0;
			for (k=n; k>=1; k--) {
				h=cc*un-un1+a[k];
				un1=un;
				un=h;
			}
			return (a[0]-un1+(c+sin(theta))*un);
		}
	}
}
