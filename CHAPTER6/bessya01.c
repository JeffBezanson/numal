#include "../real.h"


void bessya01(real_t a, real_t x, real_t *ya, real_t *ya1)
{
	if (a == 0.0) {
		void bessy01(real_t, real_t *, real_t *);
		bessy01(x,ya,ya1);
	} else {
		int n,na,rec,rev;
		real_t b,c,d,e,f,g,h,p,pi,q,r,s;
		pi=4.0*atan(1.0);
		na=floor(a+0.5);
		rec = (a >= 0.5);
		rev = (a < -0.5);
		if (rev || rec) a -= na;
		if (a == -0.5) {
			p=sqrt(2.0/pi/x);
			f=p*sin(x);
			g = -p*cos(x);
		} else if (x < 3.0) {
			real_t recipgamma(real_t, real_t *, real_t *);
			b=x/2.0;
			d = -log(b);
			e=a*d;
			c = (fabs(a) < 1.0e-8) ? 1.0/pi : a/sin(a*pi);
			s = (fabs(e) < 1.0e-8) ? 1.0 : sinh(e)/e;
			e=exp(e);
			g=recipgamma(a,&p,&q)*e;
			e=(e+1.0/e)/2.0;
			f=2.0*c*(p*e+q*s*d);
			e=a*a;
			p=g*c;
			q=1.0/g/pi;
			c=a*pi/2.0;
			r = (fabs(c) < 1.0e-8) ? 1.0 : sin(c)/c;
			r *= pi*c*r;
			c=1.0;
			d = -b*b;
			*ya = f+r*q;
			*ya1 = p;
			n=1;
			do {
				f=(f*n+p+q)/(n*n-e);
				c=c*d/n;
				p /= (n-a);
				q /= (n+a);
				g=c*(f+r*q);
				h=c*p-n*g;
				*ya += g;
				*ya1 += h;
				n++;
			} while (fabs(g/(1.0+fabs(*ya)))+fabs(h/(1.0+fabs(*ya1))) >
							1.0e-15);
			f = -(*ya);
			g = -(*ya1)/b;
		} else {
			void besspqa01(real_t,real_t,real_t *,real_t *,real_t *,real_t *);
			b=x-pi*(a+0.5)/2.0;
			c=cos(b);
			s=sin(b);
			d=sqrt(2.0/x/pi);
			besspqa01(a,x,&p,&q,&b,&h);
			f=d*(p*s+q*c);
			g=d*(h*s-b*c);
		}
		if (rev) {
			x=2.0/x;
			na = -na-1;
			for (n=0; n<=na; n++) {
				h=x*(a-n)*f-g;
				g=f;
				f=h;
			}
		} else if (rec) {
			x=2.0/x;
			for (n=1; n<=na; n++) {
				h=x*(a+n)*g-f;
				f=g;
				g=h;
			}
		}
		*ya = f;
		*ya1 = g;
	}
}
