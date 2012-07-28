#include "../real.h"


void femlagspher(real_t x[], real_t y[], int n, int nc,
			real_t (*r)(real_t), real_t (*f)(real_t), int order, real_t e[])
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int l,l1;
	real_t xl1,xl,h,a12,b1,b2,tau1,tau2,ch,tl,g,yl,pp,tau3,b3,a13,a22,
			a23,c32,c12,e1,e2,e3,e4,e5,e6,*t,*sub,*chi,*gi,xm,vl,vr,wl,
			wr,pr,rm,fm,xl2,xlxr,xr2,xlm,xrm,vlm,vrm,wlm,wrm,flm,frm,rlm,
			rrm,pl1,pl2,pl3,pr1,pr2,pr3,ql1,ql2,ql3,rlmpl1,rlmpl2,rrmpr1,
			rrmpr2,vlmql1,vlmql2,vrmqr1,vrmqr2,qr1,qr2,qr3,a,a2,a3,a4,
			b,b4,p4h,p2,p3,p4,aux1,aux2,a5,a6,a7,a8,b5,b6,b7,b8,ab4,
			a2b3,a3b2,a4b,p5,p8,p8h,aux,plm,prm;

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
			if (nc == 0)
				vl=vr=0.5;
			else if (nc == 1) {
				vl=(xl1*2.0+xl)/6.0;
				vr=(xl1+xl*2.0)/6.0;
			} else {
				xl2=xl1*xl1/12.0;
				xlxr=xl1*xl/6.0;
				xr2=xl*xl/12.0;
				vl=3.0*xl2+xlxr+xr2;
				vr=3.0*xr2+xlxr+xl2;
			}
			wl=h*vl;
			wr=h*vr;
			pr=vr/(vl+vr);
			xm=xl1+h*pr;
			fm=(*f)(xm);
			rm=(*r)(xm);
			tau1=wl*rm;
			tau2=wr*rm;
			b1=wl*fm;
			b2=wr*fm;
			a12 = -(vl+vr)/h+h*(1.0-pr)*pr*rm;
		} else {
			/* element mat vec evaluation 2 */
			if (nc == 0) {
				xlm=xl1+h*0.2113248654052;
				xrm=xl1+xl-xlm;
				vlm=vrm=0.5;
				pl1=pr3=0.45534180126148;
				pl3=pr1 = -0.12200846792815;
				pl2=pr2=1.0-pl1-pl3;
				ql1 = -2.15470053837925;
				ql3 = -0.15470053837925;
				ql2 = -ql1-ql3;
				qr1 = -ql3;
				qr3 = -ql1;
				qr2 = -ql2;
			} else if (nc == 1) {
				a=xl1;
				a2=a*a;
				a3=a*a2;
				a4=a*a3;
				b=xl;
				b2=b*b;
				b3=b*b2;
				b4=b*b3;
				p2=10.0*(a2+4.0*a*b+b2);
				p3=6.0*(a3+4.0*(a2*b+a*b2)+b3);
				p4=sqrt(6.0*(a4+10.0*(a*b3+a3*b)+28.0*a2*b2+b4));
				p4h=p4*h;
				xlm=(p3-p4h)/p2;
				xrm=(p3+p4h)/p2;
				aux1=(a+b)/4.0;
				aux2=h*(a2+7.0*a*b+b2)/6.0/p4;
				vlm=aux1-aux2;
				vrm=aux1+aux2;
			} else {
				a=xl1;
				a2=a*a;
				a3=a*a2;
				a4=a*a3;
				a5=a*a4;
				a6=a*a5;
				a7=a*a6;
				a8=a*a7;
				b=xl;
				b2=b*b;
				b3=b*b2;
				b4=b*b3;
				b5=b*b4;
				b6=b*b5;
				b7=b*b6;
				b8=b*b7;
				ab4=a*b4;
				a2b3=a2*b3;
				a3b2=a3*b2;
				a4b=a4*b;
				p4=15.0*(a4+4.0*(a3*b+a*b3)+10.0*a2*b2+b4);
				p5=10.0*(a5+4.0*(a4b+ab4)+10.0*(a3b2+a2b3)+b5);
				p8=sqrt(10.0*(a8+10.0*(a7*b+a*b7)+55.0*(a2*b6+a6*b2)+
						164.0*(a5*b3+a3*b5)+290.0*a4*b4+b8));
				aux1=(a2+a*b+b2)/6.0;
				p8h=p8*h;
				aux2=(h*(a5+7.0*(a4b+ab4)+28.0*(a3b2+a2b3)+b5))/4.8/p8;
				xlm=(p5-p8h)/p4;
				xrm=(p5+p8h)/p4;
				vlm=aux1-aux2;
				vrm=aux1+aux2;
			}
			if (nc > 0) {
				plm=(xlm-xl1)/h;
				prm=(xrm-xl1)/h;
				aux=2.0*plm-1.0;
				pl1=aux*(plm-1.0);
				pl3=aux*plm;
				pl2=1.0-pl1-pl3;
				aux=2.0*prm-1.0;
				pr1=aux*(prm-1.0);
				pr3=aux*prm;
				pr2=1.0-pr1-pr3;
				aux=4.0*plm;
				ql1=aux-3.0;
				ql3=aux-1.0;
				ql2 = -ql1-ql3;
				aux=4.0*prm;
				qr1=aux-3.0;
				qr3=aux-1.0;
				qr2 = -qr1-qr3;
			}
			wlm=h*vlm;
			wrm=h*vrm;
			vlm /= h;
			vrm /= h;
			flm=(*f)(xlm)*wlm;
			frm=wrm*(*f)(xrm);
			rlm=(*r)(xlm);
			rrm=wrm*(*r)(xrm);
			tau1=pl1*rlm+pr1*rrm;
			tau2=pl2*rlm+pr2*rrm;
			tau3=pl3*rlm+pr3*rrm;
			b1=pl1*flm+pr1*frm;
			b2=pl2*flm+pr2*frm;
			b3=pl3*flm+pr3*frm;
			vlmql1=ql1*vlm;
			vrmqr1=qr1*vrm;
			vlmql2=ql2*vlm;
			vrmqr2=qr2*vrm;
			rlmpl1=rlm*pl1;
			rrmpr1=rrm*pr1;
			rlmpl2=rlm*pl2;
			rrmpr2=rrm*pr2;
			a12=vlmql1*ql2+vrmqr1*qr2+rlmpl1*pl2+rrmpr1*pr2;
			a13=vlmql1*ql3+vrmqr1*qr3+rlmpl1*pl3+rrmpr1*pr3;
			a22=vlmql2*ql2+vrmqr2*qr2+rlmpl2*pl2+rrmpr2*pr2;
			a23=vlmql2*ql3+vrmqr2*qr3+rlmpl2*pl3+rrmpr2*pr3;
			c12 = -a12/a22;
			c32 = -a23/a22;
			a12=a13+c32*a12;
			b1 += c12*b2;
			b2=b3+c32*b2;
			tau1 += c12*tau2;
			tau2=tau3+c32*tau2;
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
				aux=((nc == 0) ? 1.0 : pow(x[0],nc))/e2;
				b1 -= e3*aux;
				tau1 -= e1*aux;
			} else if (l == n && e5 == 0.0) {
				tau2=1.0;
				b2=e6/e4;
				b1 -= a12*b2;
				tau1 -= a12;
				a12=0.0;
			} else if (l == n && e5 != 0.0) {
				aux=((nc == 0) ? 1.0 : pow(x[n],nc))/e5;
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
