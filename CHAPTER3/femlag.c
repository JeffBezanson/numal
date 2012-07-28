#include "../real.h"


void femlag(real_t x[], real_t y[], int n, real_t (*r)(real_t),
				real_t (*f)(real_t), int order, real_t e[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int l,l1;
	real_t xl1,xl,h,a12,b1,b2,tau1,tau2,ch,tl,g,yl,pp,e1,e2,e3,e4,e5,
			e6,*t,*sub,*chi,*gi,f2,r2,r1,f1,h2,r3,f3,x2,h6,h15,b3,tau3,
			c12,a13,a22,a23,r4,f4,x3,h12,h24,det,c13,c42,c43,a14,a24,
			a33,a34,b4,tau4;

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
				f2=(*f)(xl1);
				r2=(*r)(xl1);
			}
			a12 = -1.0/h;
			h2=h/2.0;
			r1=r2;
			r2=(*r)(xl);
			f1=f2;
			f2=(*f)(xl);
			b1=h2*f1;
			b2=h2*f2;
			tau1=h2*r1;
			tau2=h2*r2;
		} else if (order == 4) {
			/* element mat vec evaluation 2 */
			if (l == 1) {
				r3=(*r)(xl1);
				f3=(*f)(xl1);
			}
			x2=(xl1+xl)/2.0;
			h6=h/6.0;
			h15=h/1.5;
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
			a12 = a23 = -8.0/h/3.0;
			a13 = -a12/8.0;
			a22 = -2.0*a12+tau2;
			c12 = -a12/a22;
			a12=a13+c12*a12;
			b2 *= c12;
			b1 += b2;
			b2 += b3;
			tau2 *= c12;
			tau1 += tau2;
			tau2=tau3+tau2;
		} else {
			/* element mat vec evaluation 3 */
			if (l == 1) {
				r4=(*r)(xl1);
				f4=(*f)(xl1);
			}
			x2=xl1+0.27639320225*h;
			x3=xl-x2+xl1;
			r1=r4;
			r2=(*r)(x2);
			r3=(*r)(x3);
			r4=(*r)(xl);
			f1=f4;
			f2=(*f)(x2);
			f3=(*f)(x3);
			f4=(*f)(xl);
			h12=h/12.0;
			h24=h/2.4;
			b1=h12*f1;
			b2=h24*f2;
			b3=h24*f3;
			b4=h12*f4;
			tau1=h12*r1;
			tau2=h24*r2;
			tau3=h24*r3;
			tau4=h12*r4;
			a12 = a34 = -4.8784183052078/h;
			a13=a24=0.7117516385412/h;
			a14 = -0.16666666666667/h;
			a23=25.0*a14;
			a22 = -2.0*a23+tau2;
			a33 = -2.0*a23+tau3;
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
				tau1 -= e1/e2;
				b1 -= e3/e2;
			} else if (l == n && e5 == 0.0) {
				tau2=1.0;
				b2=e6/e4;
				b1 -= a12*b2;
				tau1 -= a12;
				a12=0.0;
			} else if (l == n && e5 != 0.0) {
				tau2 += e4/e5;
				b2 += e6/e5;
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
