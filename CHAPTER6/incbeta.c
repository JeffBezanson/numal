#include "../real.h"


real_t incbeta(real_t x, real_t p, real_t q, real_t eps)
{
	real_t numal_gamma(real_t);
	int m,n,neven,recur;
	real_t g,f,fn,fn1,fn2,gn,gn1,gn2,dn,pq;

	if (x == 0.0 || x == 1.0)
		return x;
	else {
		if (x > 0.5) {
			f=p;
			p=q;
			q=f;
			x=1.0-x;
			recur=1;
		} else
			recur=0;
		g=fn2=0.0;
		m=0;
		pq=p+q;
		f=fn1=gn1=gn2=1.0;
		neven=0;
		n=1;
		do {
			if (neven) {
				m++;
				dn=m*x*(q-m)/(p+n-1.0)/(p+n);
			} else
				dn = -x*(p+m)*(pq+m)/(p+n-1.0)/(p+n);
			g=f;
			fn=fn1+dn*fn2;
			gn=gn1+dn*gn2;
			neven=(!neven);
			f=fn/gn;
			fn2=fn1;
			fn1=fn;
			gn2=gn1;
			gn1=gn;
			n++;
		} while (fabs((f-g)/f) > eps);
		f=f*pow(x,p)*pow(1.0-x,q)*numal_gamma(p+q)/numal_gamma(p+1.0)/numal_gamma(q);
		if (recur) f=1.0-f;
		return f;
	}
}
