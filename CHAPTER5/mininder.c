#include "../real.h"


real_t mininder(real_t *x, real_t *y, real_t (*fx)(real_t),
					real_t (*dfx)(real_t), real_t (*tolx)(real_t))
{
	int sgn;
	real_t a,b,c,fa,fb,fu,dfa,dfb,dfu,e,d,tol,ba,z,p,q,s;

	if (*x <= *y) {
		a = *x;
		fa=(*fx)(*x);
		dfa=(*dfx)(*x);
		b = *x = *y;
		fb=(*fx)(*x);
		dfb=(*dfx)(*x);
	} else {
		b = *x;
		fb=(*fx)(*x);
		dfb=(*dfx)(*x);
		a = *x = *y;
		fa=(*fx)(*x);
		dfa=(*dfx)(*x);
	}
	c=(3.0-sqrt(5.0))/2.0;
	d=b-a;
	e=d*2.0;
	z=e*2.0;
	while (1) {
		ba=b-a;
		tol=(*tolx)(*x);
		if (ba < tol*3.0) break;
		if (fabs(dfa) <= fabs(dfb)) {
			*x=a;
			sgn=1;
		} else {
			*x=b;
			sgn = -1;
		}
		if (dfa <= 0.0 && dfb >= 0.0) {
			z=(fa-fb)*3.0/ba+dfa+dfb;
			s=sqrt(z*z-dfa*dfb);
			p = (sgn == 1) ? dfa-s-z : dfb+s-z;
			p *= ba;
			q=dfb-dfa+s*2.0;
			z=e;
			e=d;
			d = (fabs(p) <= fabs(q)*tol) ? tol*sgn : -p/q;
		} else
			d=ba;
		if (fabs(d) >= fabs(z*0.5) || fabs(d) > ba*0.5) {
			e=ba;
			d=c*ba*sgn;
		}
		*x += d;
		fu=(*fx)(*x);
		dfu=(*dfx)(*x);
		if (dfu >= 0.0 || (fu >= fa && dfa <= 0.0)) {
			b = *x;
			fb=fu;
			dfb=dfu;
		} else {
			a = *x;
			fa=fu;
			dfa=dfu;
		}
	}
	if (fa <= fb) {
		*x=a;
		*y=b;
		return fa;
	} else {
		*x=b;
		*y=a;
		return fb;
	}
}
