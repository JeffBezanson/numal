#include "../real.h"


void rk3(real_t *x, real_t a, real_t b, real_t *y, real_t ya,
			real_t *z, real_t za, real_t (*fxy)(real_t, real_t),
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
	first=reject=1;
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
		if (reject) {
			*x=xl;
			*y=yl;
			k0=(*fxy)(*x,*y)*h;
		} else
			k0=k5*h/hl;
		*x=xl+0.276393202250021*h;
		*y=yl+(zl*0.276393202250021+k0*0.038196601125011)*h;
		k1=(*fxy)(*x,*y)*h;
		*x=xl+0.723606797749979*h;
		*y=yl+(zl*0.723606797749979+k1*0.261803398874989)*h;
		k2=(*fxy)(*x,*y)*h;
		*x=xl+h*0.5;
		*y=yl+(zl*0.5+k0*0.046875+k1*0.079824155839840-
					k2*0.001699155839840)*h;
		k4=(*fxy)(*x,*y)*h;
		*x = (last ? b : xl+h);
		*y=yl+(zl+k0*0.309016994374947+k2*0.190983005625053)*h;
		k3=(*fxy)(*x,*y)*h;
		*y=yl+(zl+k0*0.083333333333333+k1*0.301502832395825+
					k2*0.115163834270842)*h;
		k5=(*fxy)(*x,*y)*h;
		discry=fabs((-k0*0.5+k1*1.809016994374947+
					k2*0.690983005625053-k4*2.0)*h);
		discrz=fabs((k0-k3)*2.0-(k1+k2)*10.0+k4*16.0+k5*4.0);
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
			*z=zl+(k0+k3)*0.083333333333333+
					(k1+k2)*0.416666666666667;
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
