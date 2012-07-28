#include "../real.h"


void besspqa01(real_t a, real_t x, real_t *pa, real_t *qa, real_t *pa1,
				real_t *qa1)
{
	if (a == 0.0) {
		void besspq0(real_t, real_t *, real_t *);
		void besspq1(real_t, real_t *, real_t *);
		besspq0(x,pa,qa);
		besspq1(x,pa1,qa1);
	} else {
		int n,na,rec,rev;
		real_t b,pi,p0,q0;
		pi=4.0*atan(1.0);
		rev = (a < -0.5);
		if (rev) a = -a-1.0;
		rec = (a >= 0.5);
		if (rec) {
			na=floor(a+0.5);
			a -= na;
		}
		if (a == -0.5) {
			*pa = *pa1 = 1.0;
			*qa = *qa1 = 0.0;
		} else if (x >= 3.0) {
			real_t c,d,e,f,g,p,q,r,s,temp;
			c=0.25-a*a;
			b=x+x;
			f=r=1.0;
			g = -x;
			s=0.0;
			temp=x*cos(a*pi)/pi*1.0e15;
			e=temp*temp;
			n=2;
			do {
				d=(n-1+c/n);
				p=(2*n*f+b*g-d*r)/(n+1);
				q=(2*n*g-b*f-d*s)/(n+1);
				r=f;
				f=p;
				s=g;
				g=q;
				n++;
			} while ((p*p+q*q)*n*n < e);
			e=f*f+g*g;
			p=(r*f+s*g)/e;
			q=(s*f-r*g)/e;
			f=p;
			g=q;
			n--;
			while (n > 0) {
				r=(n+1)*(2.0-p)-2.0;
				s=b+(n+1)*q;
				d=(n-1+c/n)/(r*r+s*s);
				p=d*r;
				q=d*s;
				e=f;
				f=p*(e+1.0)-g*q;
				g=q*(e+1.0)+p*g;
				n--;
			}
			f += 1.0;
			d=f*f+g*g;
			*pa = f/d;
			*qa = -g/d;
			d=a+0.5-p;
			q += x;
			*pa1 = ((*pa)*q-(*qa)*d)/x;
			*qa1 = ((*qa)*q+(*pa)*d)/x;
		} else {
			void bessjaplusn(real_t, real_t, int, real_t []);
			void bessya01(real_t, real_t, real_t *, real_t *);
			real_t c,s,chi,ya,ya1,ja[2];
			b=sqrt(pi*x/2.0);
			chi=x-pi*(a/2.0+0.25);
			c=cos(chi);
			s=sin(chi);
			bessya01(a,x,&ya,&ya1);
			bessjaplusn(a,x,1,ja);
			*pa = b*(ya*s+c*ja[0]);
			*qa = b*(c*ya-s*ja[0]);
			*pa1 = b*(s*ja[1]-c*ya1);
			*qa1 = b*(c*ja[1]+s*ya1);
		}
		if (rec) {
			x=2.0/x;
			b=(a+1.0)*x;
			for (n=1; n<=na; n++) {
				p0=(*pa)-(*qa1)*b;
				q0=(*qa)+(*pa1)*b;
				*pa = *pa1;
				*pa1 = p0;
				*qa = *qa1;
				*qa1 = q0;
				b += x;
			}
		}
		if (rev) {
			p0 = *pa1;
			*pa1 = *pa;
			*pa = p0;
			q0 = *qa1;
			*qa1 = *qa;
			*qa = q0;
		}
	}
}
