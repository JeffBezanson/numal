#include "../real.h"


void besszeros(real_t a, int n, real_t z[], int d)
{
	void besspqa01(real_t, real_t, real_t *, real_t *, real_t *, real_t *);
	int j,s;
	real_t aa,a2,b,bb,c,chi,co,mu,mu2,mu3,mu4,p,pi,pa,pa1,p0,p1,pp1,
			q,qa,qa1,q1,qq1,ro,si,t,tt,u,v,w,x,xx,x4,y,yy,fi;

	pi=4.0*atan(1.0);
	aa=a*a;
	mu=4.0*aa;
	mu2=mu*mu;
	mu3=mu*mu2;
	mu4=mu2*mu2;
	if (d < 3) {
		p=7.0*mu-31.0;
		p0=mu-1.0;
		p1=4.0*(253.0*mu2-3722.0*mu+17869.0)/15.0/p*p0;
		q1=8.0*(83.0*mu2-982.0*mu+3779.0)/5.0/p;
	} else {
		p=7.0*mu2+82.0*mu-9.0;
		p0=mu+3.0;
		p1=(4048.0*mu4+131264.0*mu3-221984.0*mu2-
			417600.0*mu+1012176.0)/60.0/p;
		q1=1.6*(83.0*mu3+2075.0*mu2-3039.0*mu+3537.0)/p;
	}
	t = (d == 1 || d == 4) ? 0.25 : 0.75;
	tt=4.0*t;
	if (d < 3) {
		pp1=5.0/48.0;
		qq1 = -5.0/36.0;
	} else {
		pp1 = -7.0/48.0;
		qq1=35.0/288.0;
	}
	y=3.0*pi/8.0;
	bb = (a >= 3.0) ? pow(a,-2.0/3.0) : 0.0;
	for (s=1; s<=n; s++) {
		if (a == 0.0 && s == 1 && d == 3) {
			x=0.0;
			j=0;
		} else {
			if (s >= 3.0*a-8.0) {
				b=(s+a/2.0-t)*pi;
				c=1.0/b/b/64.0;
				x=b-1.0/b/8.0*(p0-p1*c)/(1.0-q1*c);
			} else {
				if (s == 1)
					x = ((d == 1) ? -2.33811 : ((d == 2) ? -1.17371 :
							((d == 3) ? -1.01879 : -2.29444)));
				else {
					x=y*(4.0*s-tt);
					v=1.0/x/x;
					x = -pow(x,2.0/3.0)*(1.0+v*(pp1+qq1*v));
				}
				u=x*bb;
				yy=2.0/3.0*pow(-u,1.5);
				if (yy == 0.0)
					fi=0.0;
				else if (yy > 1.0e5)
					fi=1.570796;
				else {
					real_t r,p,pp;
					if (yy <1.0) {
						p=pow(3.0*yy,1.0/3.0);
						pp=p*p;
						p *= (1.0+pp*(-210.0+pp*(27.0-2.0*pp))/1575.0);
					} else {
						p=1.0/(yy+1.570796);
						pp=p*p;
						p=1.570796-p*(1.0+pp*(2310.0+pp*(3003.0+pp*
							(4818.0+pp*(8591.0+pp*16328.0))))/3465.0);
					}
					pp=(yy+p)*(yy+p);
					r=(p-atan(p+yy))/pp;
					fi=p-(1.0+pp)*r*(1.0+r/(p+yy));
				}
				v=fi;
				w=1.0/cos(v);
				xx=1.0-w*w;
				c=sqrt(u/xx);
				x=w*(a+c/a/u*((d < 3) ?
					-5.0/48.0/u-c*(-5.0/24.0/xx+1.0/8.0) :
					7.0/48.0/u+c*(-7.0/24.0/xx+3.0/8.0)));
			}
			j=0;
			do {
				xx=x*x;
				x4=xx*xx;
				a2=aa-xx;
				besspqa01(a,x,&pa,&qa,&pa1,&qa1);
				chi=x-pi*(a/2.0+0.25);
				si=sin(chi);
				co=cos(chi);
				ro = ((d == 1) ? (pa*co-qa*si)/(pa1*si+qa1*co) :
						((d == 2) ? (pa*si+qa*co)/(qa1*si-pa1*co) :
						((d == 3) ? a/x-(pa1*si+qa1*co)/(pa*co-qa*si) :
										a/x-(qa1*si-pa1*co)/(pa*si+qa*co))));
				j++;
				if (d < 3) {
					u=ro;
					p=(1.0-4.0*a2)/6.0/x/(2.0*a+1.0);
					q=(2.0*(xx-mu)-1.0-6.0*a)/3.0/x/(2.0*a+1.0);
				} else {
					u = -xx*ro/a2;
					v=2.0*x*a2/(aa+xx)/3.0;
					w=a2*a2*a2;
					q=v*(1.0+(mu2+32.0*mu*xx+48.0*x4)/32.0/w);
					p=v*(1.0+(-mu2+40.0*mu*xx+48.0*x4)/64.0/w);
				}
				w=u*(1.0+p*ro)/(1.0+q*ro);
				x += w;
			} while (fabs(w/x) > 1.0e-13 && j < 5);
		}
		z[s]=x;
	}
}
