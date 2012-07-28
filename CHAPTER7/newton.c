#include "../real.h"
void newton(int n, real_t x[], real_t f[])
{
	int k,i,im1;
	real_t xim1,fim1;

	im1=0;
	for (i=1; i<=n; i++) {
		fim1=f[im1];
		xim1=x[im1];
		for (k=i; k<=n; k++) f[k]=(f[k]-fim1)/(x[k]-xim1);
		im1=i;
	}
}
