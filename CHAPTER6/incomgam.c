#include "../real.h"


void incomgam(real_t x, real_t a, real_t *klgam, real_t *grgam,
				real_t gam, real_t eps)
{
	int n;
	real_t c0,c1,c2,d0,d1,d2,x2,ax,p,q,r,s,r1,r2,scf;

	s=exp(-x+a*log(x));
	scf=FLT_MAX;
	if (x <= ((a < 3.0) ? 1.0 : a)) {
		x2=x*x;
		ax=a*x;
		d0=1.0;
		p=a;
		c0=s;
		d1=(a+1.0)*(a+2.0-x);
		c1=((a+1.0)*(a+2.0)+x)*s;
		r2=c1/d1;
		n=1;
		do {
			p += 2.0;
			q=(p+1.0)*(p*(p+2.0)-ax);
			r=n*(n+a)*(p+2.0)*x2;
			c2=(q*c1+r*c0)/p;
			d2=(q*d1+r*d0)/p;
			r1=r2;
			r2=c2/d2;
			c0=c1;
			c1=c2;
			d0=d1;
			d1=d2;
			if (fabs(c1) > scf || fabs(d1) > scf) {
				c0 /= scf;
				c1 /= scf;
				d0 /= scf;
				d1 /= scf;
			}
			n++;
		} while (fabs((r2-r1)/r2) > eps);
		*klgam = r2/a;
		*grgam = gam-(*klgam);
	} else {
		c0=a*s;
		c1=(1.0+x)*c0;
		q=x+2.0-a;
		d0=x;
		d1=x*q;
		r2=c1/d1;
		n=1;
		do {
			q += 2.0;
			r=n*(n+1-a);
			c2=q*c1-r*c0;
			d2=q*d1-r*d0;
			r1=r2;
			r2=c2/d2;
			c0=c1;
			c1=c2;
			d0=d1;
			d1=d2;
			if (fabs(c1) > scf || fabs(d1) > scf) {
				c0 /= scf;
				c1 /= scf;
				d0 /= scf;
				d1 /= scf;
			}
			n++;
		} while (fabs((r2-r1)/r2) > eps);
		*grgam = r2/a;
		*klgam = gam-(*grgam);
	}
}
