#include "../real.h"


void rk2n(real_t *x, real_t a, real_t b, real_t y[], real_t ya[],
			real_t z[], real_t za[],
			real_t (*fxyzj)(int, int, real_t, real_t [], real_t []),
			real_t e[], real_t d[], int fi, int n)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int j,jj,last,first,reject,test,ta,tb;
	real_t xl,h,ind,hmin,hl,absh,fhm,discry,discrz,toly,tolz,
			mu,mu1,fhy,fhz,*yl,*zl,*k0,*k1,*k2,*k3,*k4,*k5,*ee;

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
				absh=fabs(h);
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
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj];
			z[jj]=zl[jj];
		}
		for (j=1; j<=n; j++)	k0[j]=(*fxyzj)(n,j,*x,y,z)*h;
		*x=xl+h/4.5;
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj]+(zl[jj]*18.0+k0[jj]*2.0)/81.0*h;
			z[jj]=zl[jj]+k0[jj]/4.5;
		}
		for (j=1; j<=n; j++) k1[j]=(*fxyzj)(n,j,*x,y,z)*h;
		*x=xl+h/3.0;
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj]+(zl[jj]*6.0+k0[jj])/18.0*h;
			z[jj]=zl[jj]+(k0[jj]+k1[jj]*3.0)/12.0;
		}
		for (j=1; j<=n; j++) k2[j]=(*fxyzj)(n,j,*x,y,z)*h;
		*x=xl+h*0.5;
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj]+(zl[jj]*8.0+k0[jj]+k2[jj])/16.0*h;
			z[jj]=zl[jj]+(k0[jj]+k2[jj]*3.0)/8.0;
		}
		for (j=1; j<=n; j++) k3[j]=(*fxyzj)(n,j,*x,y,z)*h;
		*x=xl+h*0.8;
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj]+(zl[jj]*100.0+k0[jj]*12.0+
						k3[jj]*28.0)/125.0*h;
			z[jj]=zl[jj]+(k0[jj]*53.0-k1[jj]*135.0+k2[jj]*126.0+
						k3[jj]*56.0)/125.0;
		}
		for (j=1; j<=n; j++) k4[j]=(*fxyzj)(n,j,*x,y,z)*h;
		*x = (last ? b : xl+h);
		for (jj=1; jj<=n; jj++) {
			y[jj]=yl[jj]+(zl[jj]*336.0+k0[jj]*21.0+k2[jj]*92.0+
						k4[jj]*55.0)/336.0*h;
			z[jj]=zl[jj]+(k0[jj]*133.0-k1[jj]*378.0+k2[jj]*276.0+
						k3[jj]*112.0+k4[jj]*25.0)/168.0;
		}
		for (j=1; j<=n; j++) k5[j]=(*fxyzj)(n,j,*x,y,z)*h;
		reject=0;
		fhm=0.0;
		for (jj=1; jj<=n; jj++) {
			discry=fabs((-k0[jj]*21.0+k2[jj]*108.0-k3[jj]*112.0+
						k4[jj]*25.0)/56.0*h);
			discrz=fabs(k0[jj]*21.0-k2[jj]*162.0+k3[jj]*224.0-
						k4[jj]*125.0+k5[jj]*42.0)/14.0;
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
					yl[jj] = y[jj];
					zl[jj] = z[jj];
				}
			} else
				h *= mu;
		} else {
			if (first) {
				first=0;
				hl=h;
				h *= mu;
			} else {
				fhm=mu*h/hl+mu-mu1;
				hl=h;
				h *= fhm;
			}
			mu1=mu;
			for (jj=1; jj<=n; jj++) {
				y[jj]=yl[jj]+(zl[jj]*56.0+k0[jj]*7.0+k2[jj]*36.0-
							k4[jj]*15.0)/56.0*hl;
				z[jj]=zl[jj]+(-k0[jj]*63.0+k1[jj]*189.0-k2[jj]*36.0-
							k3[jj]*112.0+k4[jj]*50.0)/28.0;
			}
			for (j=1; j<=n; j++) k5[j]=(*fxyzj)(n,j,*x,y,z)*hl;
			for (jj=1; jj<=n; jj++) {
				y[jj]=yl[jj]+(zl[jj]*336.0+k0[jj]*35.0+k2[jj]*108.0+
							k4[jj]*25.0)/336.0*hl;
				z[jj]=zl[jj]+(k0[jj]*35.0+k2[jj]*162.0+k4[jj]*125.0+
							k5[jj]*14.0)/336.0;
			}
			if (b == *x) break;
			xl = *x;
			for (jj=1; jj<=n; jj++) {
				yl[jj] = y[jj];
				zl[jj] = z[jj];
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
