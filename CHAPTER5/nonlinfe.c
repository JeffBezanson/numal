#include "../real.h"


void nonlinfemlagskew(real_t x[], real_t y[], int n,
				real_t (*f)(real_t, real_t, real_t),
				real_t (*fy)(real_t, real_t, real_t),
				real_t (*fz)(real_t, real_t, real_t), int nc, real_t e[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void dupvec(int, int, int, real_t [], real_t []);
	int l,l1,it;
	real_t xl1,xl,h,a12,a21,b1,b2,tau1,tau2,ch,tl,g,yl,pp,
			zl1,zl,e1,e2,e4,e5,eps,rho,*t,*super,*sub,*chi,*gi,*z,
			xm,vl,vr,wl,wr,pr,qm,rm,fm,xl12,xl1xl,xl2,zm,zaccm;

	t=allocate_real_vector(0,n-1);
	super=allocate_real_vector(0,n-1);
	sub=allocate_real_vector(0,n-1);
	chi=allocate_real_vector(0,n-1);
	gi=allocate_real_vector(0,n-1);
	z=allocate_real_vector(0,n);

	dupvec(0,n,0,z,y);
	e1=e[1];
	e2=e[2];
	e4=e[4];
	e5=e[5];
	it=1;
	do {
		l=1;
		xl=x[0];
		zl=z[0];
		while (l <= n) {
			xl1=xl;
			l1=l-1;
			xl=x[l];
			h=xl-xl1;
			zl1=zl;
			zl=z[l];
			/* element mat vec evaluation 1 */
			if (nc == 0)
				vl=vr=0.5;
			else if (nc == 1) {
				vl=(xl1*2.0+xl)/6.0;
				vr=(xl1+xl*2.0)/6.0;
			} else {
				xl12=xl1*xl1/12.0;
				xl1xl=xl1*xl/6.0;
				xl2=xl*xl/12.0;
				vl=3.0*xl12+xl1xl+xl2;
				vr=3.0*xl2+xl1xl+xl12;
			}
			wl=h*vl;
			wr=h*vr;
			pr=vr/(vl+vr);
			xm=xl1+h*pr;
			zm=pr*zl+(1.0-pr)*zl1;
			zaccm=(zl-zl1)/h;
			qm=(*fz)(xm,zm,zaccm);
			rm=(*fy)(xm,zm,zaccm);
			fm=(*f)(xm,zm,zaccm);
			tau1=wl*rm;
			tau2=wr*rm;
			b1=wl*fm-zaccm*(vl+vr);
			b2=wr*fm+zaccm*(vl+vr);
			a12 = -(vl+vr)/h+vl*qm+(1.0-pr)*pr*rm*(wl+wr);
			a21 = -(vl+vr)/h-vr*qm+(1.0-pr)*pr*rm*(wl+wr);
         if (l == 1 || l == n) {
				/* boundary conditions */
				if (l == 1 && e2 == 0.0) {
					tau1=1.0;
					b1=a12=0.0;
				} else if (l == 1 && e2 != 0.0) {
					tau1 -= e1/e2;
				} else if (l == n && e5 == 0.0) {
					tau2=1.0;
					b2=a21=0.0;
				} else if (l == n && e5 != 0.0) {
					tau2 += e4/e5;
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
		eps=0.0;
		rho=1.0;
		for (l=0; l<=n; l++) {
			rho += fabs(z[l]);
			eps += fabs(y[l]);
			z[l] -= y[l];
		}
		rho *= FLT_EPSILON;
		it++;
	} while (eps > rho);
	dupvec(0,n,0,y,z);
   free_real_vector(t,0);
	free_real_vector(super,0);
	free_real_vector(sub,0);
	free_real_vector(chi,0);
	free_real_vector(gi,0);
	free_real_vector(z,0);
}
