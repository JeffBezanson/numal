#include "../real.h"
void intchs(int n, real_t a[], real_t b[])
{
	int i;
	real_t h,l,dum;

	if (n == 0) {
		b[1]=a[0];
		return;
	}
	if (n == 1) {
		b[2]=a[1]/4.0;
		b[1]=a[0];
		return;
	}
	h=a[n];
	dum=a[n-1];
	b[n+1]=h/((n+1)*2);
	b[n]=dum/(n*2);
	for (i=n-1; i>=2; i--) {
		l=a[i-1];
		b[i]=(l-h)/(2*i);
		h=dum;
		dum=l;
	}
	b[1]=a[0]-h/2.0;
}
