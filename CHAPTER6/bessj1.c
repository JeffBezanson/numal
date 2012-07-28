#include "../real.h"


real_t bessj1(real_t x)
{
	if (x == 0.0) return 1.0;
	if (fabs(x) < 8.0) {
		int i;
		real_t z,z2,b0,b1,b2;
		static real_t ar[15]={-0.19554e-15, 0.1138572e-13,
			-0.57774042e-12, 0.2528123664e-10, -0.94242129816e-9,
			0.2949707007278e-7, -0.76175878054003e-6,
			0.158870192399321e-4, -0.260444389348581e-3,
			0.324027018268386e-2, -0.291755248061542e-1,
			0.177709117239728e0, -0.661443934134543e0,
			0.128799409885768e1, -0.119180116054122e1};
		x /= 8.0;
		z=2.0*x*x-1.0;
		z2=z+z;
		b1=b2=0.0;
		for (i=0; i<=14; i++) {
			b0=z2*b1-b2+ar[i];
			b2=b1;
			b1=b0;
		}
		return x*(z*b1-b2+0.648358770605265);
	} else {
		void besspq1(real_t, real_t *, real_t *);
		int sgnx;
		real_t c,cosx,sinx,p1,q1;
		sgnx = (x > 0.0) ? 1 : -1;
		x=fabs(x);
		c=0.797884560802865/sqrt(x);
		cosx=cos(x-0.706858347057703e1);
		sinx=sin(x-0.706858347057703e1);
		besspq1(x,&p1,&q1);
		return sgnx*c*(p1*sinx+q1*cosx);
	}
}
