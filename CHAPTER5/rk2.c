#include "../real.h"


void rk2(real_t *x, real_t a, real_t b, real_t *y, real_t ya,
			real_t *z, real_t za, real_t (*fxyz)(real_t, real_t, real_t),
			real_t e[], real_t d[], int fi)
{
	int last,first,reject,test,ta,tb;
	real_t e1,e2,e3,e4,xl,yl,zl,h,ind,hmin,hl,absh,k0,k1,k2,k3,k4,
			k5,discry,discrz,toly,tolz,mu,mu1,fhy,fhz;

	if (fi) {
		d[3]=a;
		d[4]=ya;
		d[5]=za;
	}
	d[1]=0.0;
	xl=d[3];
	yl=d[4];
	zl=d[5];
	if (fi) d[2]=b-d[3];
	absh=h=fabs(d[2]);
	if (b-xl < 0.0) h = -h;
	ind=fabs(b-xl);
	hmin=ind*e[1]+e[2];
	hl=ind*e[3]+e[4];
	if (hl < hmin) hmin=hl;
	e1=e[1]/ind;
	e2=e[2]/ind;
	e3=e[3]/ind;
	e4=e[4]/ind;
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
		*z=zl;
		k0=(*fxyz)(*x,*y,*z)*h;
		*x=xl+h/4.5;
		*y=yl+(zl*18.0+k0*2.0)/81.0*h;
		*z=zl+k0/4.5;
		k1=(*fxyz)(*x,*y,*z)*h;
		*x=xl+h/3.0;
		*y=yl+(zl*6.0+k0)/18.0*h;
		*z=zl+(k0+k1*3.0)/12.0;
		k2=(*fxyz)(*x,*y,*z)*h;
		*x=xl+h*0.5;
		*y=yl+(zl*8.0+k0+k2)/16.0*h;
		*z=zl+(k0+k2*3.0)/8.0;
		k3=(*fxyz)(*x,*y,*z)*h;
		*x=xl+h*0.8;
		*y=yl+(zl*100.0+k0*12.0+k3*28.0)/125.0*h;
		*z=zl+(k0*53.0-k1*135.0+k2*126.0+k3*56.0)/125.0;
		k4=(*fxyz)(*x,*y,*z)*h;
		*x = (last ? b : xl+h);
		*y=yl+(zl*336.0+k0*21.0+k2*92.0+k4*55.0)/336.0*h;
		*z=zl+(k0*133.0-k1*378.0+k2*276.0+k3*112.0+k4*25.0)/168.0;
		k5=(*fxyz)(*x,*y,*z)*h;
		discry=fabs((-k0*21.0+k2*108.0-k3*112.0+k4*25.0)/56.0*h);
		discrz=fabs(k0*21.0-k2*162.0+k3*224.0-k4*125.0+k5*42.0)/14.0;
		toly=absh*(fabs(zl)*e1+e2);
		tolz=fabs(k0)*e3+absh*e4;
		reject=(discry > toly || discrz > tolz);
		fhy=discry/toly;
		fhz=discrz/tolz;
		if (fhz > fhy) fhy=fhz;
		mu=1.0/(1.0+fhy)+0.45;
		if (reject) {
			if (absh <= hmin) {
				d[1] += 1.0;
				*y=yl;
				*z=zl;
				first=1;
				if (b == *x) break;
				xl = *x;
				yl = *y;
				zl = *z;
			} else
				h *= mu;
		} else {
			if (first) {
				first=0;
				hl=h;
				h *= mu;
			} else {
				fhy=mu*h/hl+mu-mu1;
				hl=h;
				h *= fhy;
			}
			mu1=mu;
			*y=yl+(zl*56.0+k0*7.0+k2*36.0-k4*15.0)/56.0*hl;
			*z=zl+(-k0*63.0+k1*189.0-k2*36.0-k3*112.0+k4*50.0)/28.0;
			k5=(*fxyz)(*x,*y,*z)*hl;
			*y=yl+(zl*336.0+k0*35.0+k2*108.0+k4*25.0)/336.0*hl;
			*z=zl+(k0*35.0+k2*162.0+k4*125.0+k5*14.0)/336.0;
			if (b == *x) break;
			xl = *x;
			yl = *y;
			zl = *z;
		}
	}
	if (!last) d[2]=h;
	d[3] = *x;
	d[4] = *y;
	d[5] = *z;
}
