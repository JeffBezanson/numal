#include "../real.h"


void besska01(real_t a, real_t x, real_t *ka, real_t *ka1)
{
	if (a == 0.0) {
		void bessk01(real_t, real_t *, real_t *);
		bessk01(x,ka,ka1);
	} else {
		int n,na,rec,rev;
		real_t f,g,h,pi;
		pi=4.0*atan(1.0);
		rev = (a < -0.5);
		if (rev) a = -a-1.0;
		rec = (a >= 0.5);
		if (rec) {
			na=floor(a+0.5);
			a -= na;
		}
		if (a == 0.5)
			f=g=sqrt(pi/x/2.0)*exp(-x);
		else if (x < 1.0) {
			real_t recipgamma(real_t, real_t *, real_t *);
			real_t a1,b,c,d,e,p,q,s;
			b=x/2.0;
			d = -log(b);
			e=a*d;
			c=a*pi;
			c = (fabs(c) < 1.0e-15) ? 1.0 : c/sin(c);
			s = (fabs(e) < 1.0e-15) ? 1.0 : sinh(e)/e;
			e=exp(e);
			a1=(e+1.0/e)/2.0;
			g=recipgamma(a,&p,&q)*e;
			*ka = f = c*(p*a1+q*s*d);
			e=a*a;
			p=0.5*g*c;
			q=0.5/g;
			c=1.0;
			d=b*b;
			*ka1 = p;
			n=1;
			do {
				f=(f*n+p+q)/(n*n-e);
				c=c*d/n;
				p /= (n-a);
				q /= (n+a);
				g=c*(p-n*f);
				h=c*f;
				*ka += h;
				*ka1 += g;
				n++;
			} while (h/(*ka)+fabs(g)/(*ka1) > 1.0e-15);
			f=(*ka);
			g=(*ka1)/b;
		} else {
			void nonexpbesska01(real_t, real_t, real_t *, real_t *);
			real_t expon;
			expon=exp(-x);
			nonexpbesska01(a,x,ka,ka1);
			f=expon*(*ka);
			g=expon*(*ka1);
		}
		if (rec) {
			x=2.0/x;
			for (n=1; n<=na; n++) {
				h=f+(a+n)*x*g;
				f=g;
				g=h;
			}
		}
		if (rev) {
			*ka1 = f;
			*ka = g;
		} else {
			*ka = f;
			*ka1 = g;
		}
	}
}
