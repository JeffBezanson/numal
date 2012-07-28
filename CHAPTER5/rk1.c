#include "../real.h"


void rk1(real_t *x, real_t a, real_t b, real_t *y, real_t ya,
			real_t (*fxy)(real_t, real_t), real_t e[], real_t d[], int fi)
{
	int last,first,reject,test,ta,tb;
	real_t e1,e2,xl,yl,h,ind,hmin,absh,k0,k1,k2,k3,k4,k5,discr,tol,mu,
			mu1,fh,hl;

	if (fi) {
		d[3]=a;
		d[4]=ya;
	}
	d[1]=0.0;
	xl=d[3];
	yl=d[4];
	if (fi) d[2]=b-d[3];
	absh=h=fabs(d[2]);
	if (b-xl < 0.0) h = -h;
	ind=fabs(b-xl);
	hmin=ind*e[1]+e[2];
	e1=e[1]/ind;
	e2=e[2]/ind;
	first=1;
	test=1;
	if (fi) {
		last=1;
		test=0;
	}
	while (1) {
		if (test) {
			absh=fabs(h);
			if (absh < hmin) {
				h = (h > 0.0) ? hmin : -hmin;
				absh=hmin;
			}
			ta=(h >= b-xl);
			tb=(h >= 0.0);
			if ((ta && tb) || (!(ta || tb))) {
				d[2]=h;
				last=1;
				h=b-xl;
				absh=fabs(h);
			} else
				last=0;
		}
		test=1;
		*x=xl;
		*y=yl;
		k0=(*fxy)(*x,*y)*h;
		*x=xl+h/4.5;
		*y=yl+k0/4.5;
		k1=(*fxy)(*x,*y)*h;
		*x=xl+h/3.0;
		*y=yl+(k0+k1*3.0)/12.0;
		k2=(*fxy)(*x,*y)*h;
		*x=xl+h*0.5;
		*y=yl+(k0+k2*3.0)/8.0;
		k3=(*fxy)(*x,*y)*h;
		*x=xl+h*0.8;
		*y=yl+(k0*53.0-k1*135.0+k2*126.0+k3*56.0)/125.0;
		k4=(*fxy)(*x,*y)*h;
		*x = (last ? b : xl+h);
		*y=yl+(k0*133.0-k1*378.0+k2*276.0+k3*112.0+k4*25.0)/168.0;
		k5=(*fxy)(*x,*y)*h;
		discr=fabs(k0*21.0-k2*162.0+k3*224.0-k4*125.0+k5*42.0)/14.0;
		tol=fabs(k0)*e1+absh*e2;
		reject = discr > tol;
		mu=tol/(tol+discr)+0.45;
		if (reject) {
			if (absh <= hmin) {
				d[1] += 1.0;
				*y=yl;
				first=1;
				if (b == *x) break;
				xl = *x;
				yl = *y;
			} else
				h *= mu;
		} else {
			if (first) {
				first=0;
				hl=h;
				h *= mu;
			} else {
				fh=mu*h/hl+mu-mu1;
				hl=h;
				h *= fh;
			}
			mu1=mu;
			*y=yl+(-k0*63.0+k1*189.0-k2*36.0-k3*112.0+k4*50.0)/28.0;
			k5=(*fxy)(*x,*y)*hl;
			*y=yl+(k0*35.0+k2*162.0+k4*125.0+k5*14.0)/336.0;
			if (b == *x) break;
			xl = *x;
			yl = *y;
		}
	}
	if (!last) d[2]=h;
	d[3] = *x;
	d[4] = *y;
}
