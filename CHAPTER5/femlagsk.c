#include "../real.h"


void femlagskew(real_t x[], real_t y[], int n, real_t (*q)(real_t),
			real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int l,l1;
	real_t xl1,xl,h,a12,a21,b1,b2,tau1,tau2,ch,tl,g,yl,pp,e1,e2,e3,e4,
			e5,e6,*t,*super,*sub,*chi,*gi,q2,r2,f2,q1,r1,f1,h2,s12,q3,
			r3,f3,s13,s22,x2,h6,h15,c12,c32,a13,a31,a22,a23,a32,b3,
			tau3,q4,r4,f4,s14,s23,x3,h12,h24,det,c13,c42,c43,a14,a24,
			a33,a34,a41,a42,a43,b4,tau4;

	t=allocate_real_vector(0,n-1);
	super=allocate_real_vector(0,n-1);
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
		xl1=xl;
		l1=l-1;
		xl=x[l];
		h=xl-xl1;
		if (order == 2) {
			/* element mat vec evaluation 1 */
			if (l == 1) {
				q2=(*q)(xl1);
				r2=(*r)(xl1);
				f2=(*f)(xl1);
			}
			h2=h/2.0;
			s12 = -1.0/h;
			q1=q2;
			q2=(*q)(xl);
			r1=r2;
			r2=(*r)(xl);
			f1=f2;
			f2=(*f)(xl);
			b1=h2*f1;
			b2=h2*f2;
			tau1=h2*r1;
			tau2=h2*r2;
			a12=s12+q1/2.0;
			a21=s12-q2/2.0;
		} else if (order == 4) {
			/* element mat vec evaluation 2 */
			if (l == 1) {
				q3=(*q)(xl1);
				r3=(*r)(xl1);
				f3=(*f)(xl1);
			}
			x2=(xl1+xl)/2.0;
			h6=h/6.0;
			h15=h/1.5;
			q1=q3;
			q2=(*q)(x2);
			q3=(*q)(xl);
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
			s12 = -1.0/h/0.375;
			s13 = -s12/8.0;
			s22 = -2.0*s12;
			a12=s12+q1/1.5;
			a13=s13-q1/6.0;
			a21=s12-q2/1.5;
			a23=s12+q2/1.5;
			a22=s22+tau2;
			a31=s13+q3/6.0;
			a32=s12-q3/1.5;
			c12 = -a12/a22;
			c32 = -a32/a22;
			a12=a13+c12*a23;
			a21=a31+c32*a21;
			b1 += c12*b2;
			b2=b3+c32*b2;
			tau1 += c12*tau2;
			tau2=tau3+c32*tau2;
		} else {
			/* element mat vec evaluation 3 */
			if (l == 1) {
				q4=(*q)(xl1);
				r4=(*r)(xl1);
				f4=(*f)(xl1);
			}
			x2=xl1+0.27639320225*h;
			x3=xl-x2+xl1;
			h12=h/12.0;
			h24=h/2.4;
			q1=q4;
			q2=(*q)(x2);
			q3=(*q)(x3);
			q4=(*q)(xl);
			r1=r4;
			r2=(*r)(x2);
			r3=(*r)(x3);
			r4=(*r)(xl);
			f1=f4;
			f2=(*f)(x2);
			f3=(*f)(x3);
			f4=(*f)(xl);
			s12 = -4.8784183052080/h;
			s13=0.7117516385414/h;
			s14 = -0.16666666666667/h;
			s23=25.0*s14;
			s22 = -2.0*s23;
			b1=h12*f1;
			b2=h24*f2;
			b3=h24*f3;
			b4=h12*f4;
			tau1=h12*r1;
			tau2=h24*r2;
			tau3=h24*r3;
			tau4=h12*r4;
			a12=s12+0.67418082864578*q1;
			a13=s13-0.25751416197912*q1;
			a14=s14+q1/12.0;
			a21=s12-0.67418082864578*q2;
			a22=s22+tau2;
			a23=s23+0.93169499062490*q2;
			a24=s13-0.25751416197912*q2;
			a31=s13+0.25751416197912*q3;
			a32=s23-0.93169499062490*q3;
			a33=s22+tau3;
			a34=s12+0.67418082864578*q3;
			a41=s14-q4/12.0;
			a42=s13+0.25751416197912*q4;
			a43=s12-0.67418082864578*q4;
			det=a22*a33-a23*a32;
			c12=(a13*a32-a12*a33)/det;
			c13=(a12*a23-a13*a22)/det;
			c42=(a32*a43-a42*a33)/det;
			c43=(a42*a23-a43*a22)/det;
			tau1 += c12*tau2+c13*tau3;
			tau2=tau4+c42*tau2+c43*tau3;
			a12=a14+c12*a24+c13*a34;
			a21=a41+c42*a21+c43*a31;
			b1 += c12*b2+c13*b3;
			b2=b4+c42*b2+c43*b3;
		}
		if (l == 1 || l == n) {
			/* boundary conditions */
			if (l == 1 && e2 == 0.0) {
				tau1=1.0;
				b1=e3/e1;
				a12=0.0;
			} else if (l == 1 && e2 != 0.0) {
				tau1 -= e1/e2;
				b1 -= e3/e2;
			} else if (l == n && e5 == 0.0) {
				tau2=1.0;
				a21=0.0;
				b2=e6/e4;
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
			sub[0]=a21;
			super[0]=a12;
			pp=a21/(ch-a12);
			ch=tau2-ch*pp;
			g=b2-g*pp;
			tl=tau2;
			yl=b2;
		} else {
			chi[l1] = ch += tau1;
			gi[l1] = g += b1;
			sub[l1]=a21;
			super[l1]=a12;
			pp=a21/(ch-a12);
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
		pp=super[l]/(ch-sub[l]);
		tl=t[l];
		ch=tl-ch*pp;
		yl=y[l];
		g=yl-g*pp;
		y[l]=(gi[l]+g-yl)/(chi[l]+ch-tl);
		l--;
	}
	free_real_vector(t,0);
	free_real_vector(super,0);
	free_real_vector(sub,0);
	free_real_vector(chi,0);
	free_real_vector(gi,0);
}
