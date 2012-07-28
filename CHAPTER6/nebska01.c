#include "../real.h"


void nonexpbesska01(real_t a, real_t x, real_t *ka, real_t *ka1)
{
	if (a == 0.0) {
		void nonexpbessk01(real_t, real_t *, real_t *);
		nonexpbessk01(x,ka,ka1);
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
		if (a == -0.5)
			f=g=sqrt(pi/x/2.0);
		else if (x < 1.0) {
			void besska01(real_t, real_t, real_t *, real_t *);
			real_t expon;
			expon=exp(x);
			besska01(a,x,ka,ka1);
			f=expon*(*ka);
			g=expon*(*ka1);
		} else {
			real_t b,c,e,p,q;
			c=0.25-a*a;
			b=x+x;
			g=1.0;
			f=0.0;
			e=cos(a*pi)/pi*x*1.0e15;
			n=1;
			do {
				h=(2.0*(n+x)*g-(n-1+c/n)*f)/(n+1);
				f=g;
				g=h;
				n++;
			} while (h*n < e);
			p=q=f/g;
			e=b-2.0;
			do {
				p=(n-1+c/n)/(e+(n+1)*(2.0-p));
				q=p*(1.0+q);
				n--;
			} while (n > 0);
			f=sqrt(pi/b)/(1.0+q);
			g=f*(a+x+0.5-p)/x;
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
