#include "../real.h"


void rk3n(real_t *x, real_t a, real_t b, real_t y[], real_t ya[],
			real_t z[], real_t za[],
			real_t (*fxyj)(int, int, real_t, real_t[]),
			real_t e[], real_t d[], int fi, int n)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int j,jj,last,first,reject,test,ta,tb;
	real_t xl,h,hmin,ind,hl,absh,fhm,discry,discrz,toly,tolz,mu,
			mu1,fhy,fhz,*yl,*zl,*k0,*k1,*k2,*k3,*k4,*k5,*ee;

	yl=allocate_real_vector(1,n);
	zl=allocate_real_vector(1,n);
	k0=allocate_real_vector(1,n);
	k1=allocate_real_vector(1,n);
	k2=allocate_real_vector(1,n);
	k3=allocate_real_vector(1,n);
	k4=allocate_real_vector(1,n);
	k5=allocate_real_vector(1,n);
	ee=allocate_real_vector(1,4*n);

	if (fi) {
		d[3]=a;
		for (jj=1; jj<=n; jj++) {
			d[jj+3]=ya[jj];
			d[n+jj+3]=za[jj];
		}
	}
	d[1]=0.0;
	xl=d[3];
	for (jj=1; jj<=n; jj++) {
		yl[jj]=d[jj+3];
		zl[jj]=d[n+jj+3];
	}
	if (fi) d[2]=b-d[3];
	absh=h=fabs(d[2]);
	if (b-xl < 0.0) h = -h;
	ind=fabs(b-xl);
	hmin=ind*e[1]+e[2];
	for (jj=2; jj<=2*n; jj++) {
		hl=ind*e[2*jj-1]+e[2*jj];
		if (hl < hmin) hmin=hl;
	}
	for (jj=1; jj<=4*n; jj++) ee[jj]=e[jj]/ind;
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
			for (jj=1; jj<=n; jj++) y[jj]=yl[jj];
			for (j=1; j<=n; j++) k0[j]=(*fxyj)(n,j,*x,y)*h;
		} else {
			fhy=h/hl;
			for (jj=1; jj<=n; jj++) k0[jj]=k5[jj]*fhy;
		}
		*x=xl+0.276393202250021*h;
		for (jj=1; jj<=n; jj++)
			y[jj]=yl[jj]+(zl[jj]*0.276393202250021+
							k0[jj]*0.038196601125011)*h;
		for (j=1; j<=n; j++) k1[j]=(*fxyj)(n,j,*x,y)*h;
		*x=xl+0.723606797749979*h;
		for (jj=1; jj<=n; jj++)
			y[jj]=yl[jj]+(zl[jj]*0.723606797749979+
							k1[jj]*0.261803398874989)*h;
		for (j=1; j<=n; j++) k2[j]=(*fxyj)(n,j,*x,y)*h;
		*x=xl+h*0.5;
		for (jj=1; jj<=n; jj++)
			y[jj]=yl[jj]+(zl[jj]*0.5+k0[jj]*0.046875+k1[jj]*
					0.079824155839840-k2[jj]*0.001699155839840)*h;
		for (j=1; j<=n; j++) k4[j]=(*fxyj)(n,j,*x,y)*h;
		*x = (last ? b : xl+h);
		for (jj=1; jj<=n; jj++)
			y[jj]=yl[jj]+(zl[jj]+k0[jj]*0.309016994374947+
							k2[jj]*0.190983005625053)*h;
		for (j=1; j<=n; j++) k3[j]=(*fxyj)(n,j,*x,y)*h;
		for (jj=1; jj<=n; jj++)
			y[jj]=yl[jj]+(zl[jj]+k0[jj]*0.083333333333333+k1[jj]*
					0.301502832395825+k2[jj]*0.115163834270842)*h;
		for (j=1; j<=n; j++) k5[j]=(*fxyj)(n,j,*x,y)*h;
		reject=0;
		fhm=0.0;
		for (jj=1; jj<=n; jj++) {
			discry=fabs((-k0[jj]*0.5+k1[jj]*1.809016994374947+
						k2[jj]*0.690983005625053-k4[jj]*2.0)*h);
			discrz=fabs((k0[jj]-k3[jj])*2.0-(k1[jj]+k2[jj])*10.0+
						k4[jj]*16.0+k5[jj]*4.0);
			toly=absh*(fabs(zl[jj])*ee[2*jj-1]+ee[2*jj]);
			tolz=fabs(k0[jj])*ee[2*(jj+n)-1]+absh*ee[2*(jj+n)];
			reject=((discry > toly) || (discrz > tolz) || reject);
			fhy=discry/toly;
			fhz=discrz/tolz;
			if (fhz > fhy) fhy=fhz;
			if (fhy > fhm) fhm=fhy;
		}
		mu=1.0/(1.0+fhm)+0.45;
		if (reject) {
			if (absh <= hmin) {
				d[1] += 1.0;
				for (jj=1; jj<=n; jj++) {
					y[jj]=yl[jj];
					z[jj]=zl[jj];
				}
				first=1;
				if (b == *x) break;
				xl = *x;
				for (jj=1; jj<=n; jj++) {
					yl[jj]=y[jj];
					zl[jj]=z[jj];
				}
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
			for (jj=1; jj<=n; jj++)
				z[jj]=zl[jj]+(k0[jj]+k3[jj])*0.083333333333333+
								(k1[jj]+k2[jj])*0.416666666666667;
			if (b == *x) break;
			xl = *x;
			for (jj=1; jj<=n; jj++) {
				yl[jj]=y[jj];
				zl[jj]=z[jj];
			}
		}
	}
	if (!last) d[2]=h;
	d[3] = *x;
	for (jj=1; jj<=n; jj++) {
		d[jj+3]=y[jj];
		d[n+jj+3]=z[jj];
	}
	free_real_vector(yl,1);
	free_real_vector(zl,1);
	free_real_vector(k0,1);
	free_real_vector(k1,1);
	free_real_vector(k2,1);
	free_real_vector(k3,1);
	free_real_vector(k4,1);
	free_real_vector(k5,1);
	free_real_vector(ee,1);
}
