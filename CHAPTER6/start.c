#include "../real.h"


int start(real_t x, int n, int t)
{
	int s;
	real_t p,q,r,y;

	s=2*t-1;
	p=36.0/x-t;
	r=n/x;
	if (r > 1.0 || t == 1) {
		q=sqrt(r*r+s);
		r=r*log(q+r)-q;
	} else
		r=0.0;
	q=18.0/x+r;
	r = (p > q) ? p : q;
	p=sqrt(2.0*(t+r));
	p=x*((1.0+r)+p)/(1.0+p);
	y=0.0;
	q=y;
	do {
		y=p;
		p /= x;
		q=sqrt(p*p+s);
		p=x*(r+q)/log(p+q);
		q=y;
	} while (p > q || p < q-1.0);
	return ((t == 1) ? floor(p+1.0) : -floor(-p/2.0)*2);
}
