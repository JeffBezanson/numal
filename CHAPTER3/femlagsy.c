#include "../real.h"


void femlagsym(real_t x[], real_t y[], int n, real_t (*p)(real_t),
			real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int l,l1;
	real_t xl1,xl,h,a12,b1,b2,tau1,tau2,ch,tl,g,yl,pp,p1,p2,p3,p4,
			r1,r2,r3,r4,f1,f2,f3,f4,e1,e2,e3,e4,e5,e6,*t,*sub,*chi,*gi,
			h2,x2,h6,h15,b3,tau3,c12,c32,a13,a22,a23,x3,h12,h24,det,
			c13,c42,c43,a14,a24,a33,a34,b4,tau4,aux;

	t=allocate_real_vector(0,n-1);
	sub=allocate_real_vector(0,n-1);
	chi=allocate_real_vector(0,n-1);
	gi=allocate_real_vector(0,n-1);

	l=1;
	xl=x[0];
	e1=e[1];
	e2=e[2];
	e3=e[3];
	e4=e[4];
	e5=e[5];
	e6=e[6];
	while (l <= n) {
		l1=l-1;
		xl1=xl;
		xl=x[l];
		h=xl-xl1;
		if (order == 2) {
			/* element mat vec evaluation 1 */
			if (l == 1) {
				p2=(*p)(xl1);
				r2=(*r)(xl1);
				f2=(*f)(xl1);
			}
			p1=p2;
			p2=(*p)(xl);
			r1=r2;
			r2=(*r)(xl);
			f1=f2;
			f2=(*f)(xl);
			h2=h/2.0;
			b1=h2*f1;
			b2=h2*f2;
			tau1=h2*r1;
			tau2=h2*r2;
			a12 = -0.5*(p1+p2)/h;
		} else if (order == 4) {
			/* element mat vec evaluation 2 */
			if (l == 1) {
				p3=(*p)(xl1);
				r3=(*r)(xl1);
				f3=(*f)(xl1);
			}
			x2=(xl1+xl)/2.0;
			h6=h/6.0;
			h15=h/1.5;
			p1=p3;
			p2=(*p)(x2);
			p3=(*p)(xl);
			r1=r3;
			r2=(*r)(x2);
			r3=(*r)(xl);
			f1=f3;
			f2=(*f)(x2);
			f3=(*f)(xl);
			b1=h6*f1;
			b2=h15*f2;
			b3=h6*f3;
			tau1=h6*r1;
			tau2=h15*r2;
			tau3=h6*r3;
			a12 = -(2.0*p1+p3/1.5)/h;
			a13=(0.5*(p1+p3)-p2/1.5)/h;
			a22=(p1+p3)/h/0.375+tau2;
			a23 = -(p1/3.0+p3)*2.0/h;
			c12 = -a12/a22;
			c32 = -a23/a22;
			a12=a13+c32*a12;
			b1 += c12*b2;
			b2=b3+c32*b2;
			tau1 += c12*tau2;
			tau2=tau3+c32*tau2;
		} else {
			/* element mat vec evaluation 3 */
			if (l == 1) {
				p4=(*p)(xl1);
				r4=(*r)(xl1);
				f4=(*f)(xl1);
			}
			x2=xl1+0.27639320225*h;
			x3=xl-x2+xl1;
			h12=h/12.0;
			h24=h/2.4;
			p1=p4;
			p2=(*p)(x2);
			p3=(*p)(x3);
			p4=(*p)(xl);
			r1=r4;
			r2=(*r)(x2);
			r3=(*r)(x3);
			r4=(*r)(xl);
			f1=f4;
			f2=(*f)(x2);
			f3=(*f)(x3);
			f4=(*f)(xl);
			b1=h12*f1;
			b2=h24*f2;
			b3=h24*f3;
			b4=h12*f4;
			tau1=h12*r1;
			tau2=h24*r2;
			tau3=h24*r3;
			tau4=h12*r4;
			a12 = -(4.04508497187450*p1+0.57581917135425*p3+
						0.25751416197911*p4)/h;
			a13=(1.5450849718747*p1-1.5075141619791*p2+
						0.6741808286458*p4)/h;
			a14=((p2+p3)/2.4-(p1+p4)/2.0)/h;
			a22=(5.454237476562*p1+p3/0.48+0.79576252343762*p4)/h+tau2;
			a23 = -(p1+p4)/(h*0.48);
			a24=(0.67418082864575*p1-1.50751416197910*p3+
					1.54508497187470*p4)/h;
			a33=(0.7957625234376*p1+p2/0.48+5.454237476562*p4)/h+tau3;
			a34 = -(0.25751416197911*p1+0.57581917135418*p2+
					4.0450849718747*p4)/h;
			det=a22*a33-a23*a23;
			c12=(a13*a23-a12*a33)/det;
			c13=(a12*a23-a13*a22)/det;
			c42=(a23*a34-a24*a33)/det;
			c43=(a24*a23-a34*a22)/det;
			tau1 += c12*tau2+c13*tau3;
			tau2=tau4+c42*tau2+c43*tau3;
			a12=a14+c42*a12+c43*a13;
			b1 += c12*b2+c13*b3;
			b2=b4+c42*b2+c43*b3;
		}
		if (l == 1 || l == n) {
			/* boundary conditions */
			if (l == 1 && e2 == 0.0) {
				tau1=1.0;
				b1=e3/e1;
				b2 -= a12*b1;
				tau2 -= a12;
				a12=0.0;
			} else if (l == 1 && e2 != 0.0) {
				aux=p1/e2;
				tau1 -= aux*e1;
				b1 -= e3*aux;
			} else if (l == n && e5 == 0.0) {
				tau2=1.0;
				b2=e6/e4;
				b1 -= a12*b2;
				tau1 -= a12;
				a12=0.0;
			} else if (l == n && e5 != 0.0) {
				aux=p2/e5;
				tau2 += aux*e4;
				b2 += aux*e6;
			}
		}
		/* forward babushka */
		if (l == 1) {
			chi[0]=ch=tl=tau1;
			t[0]=tl;
			gi[0]=g=yl=b1;
			y[0]=yl;
			sub[0]=a12;
			pp=a12/(ch-a12);
			ch=tau2-ch*pp;
			g=b2-g*pp;
			tl=tau2;
			yl=b2;
		} else {
			chi[l1] = ch += tau1;
			gi[l1] = g += b1;
			sub[l1]=a12;
			pp=a12/(ch-a12);
			ch=tau2-ch*pp;
			g=b2-g*pp;
			t[l1]=tl+tau1;
			tl=tau2;
			y[l1]=yl+b1;
			yl=b2;
		}
		l++;
	}
	/* backward babushka */
	pp=yl;
	y[n]=g/ch;
	g=pp;
	ch=tl;
	l=n-1;
	while (l >= 0) {
		pp=sub[l];
		pp /= (ch-pp);
		tl=t[l];
		ch=tl-ch*pp;
		yl=y[l];
		g=yl-g*pp;
		y[l]=(gi[l]+g-yl)/(chi[l]+ch-tl);
		l--;
	}
	free_real_vector(t,0);
	free_real_vector(sub,0);
	free_real_vector(chi,0);
	free_real_vector(gi,0);
}
