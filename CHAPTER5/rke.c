#include "../real.h"


void rke(real_t *x, real_t *xe, int n, real_t y[],
			void (*der)(int, real_t, real_t[]), real_t data[], int fi,
			void (*out)(int, real_t, real_t, real_t [], real_t []))
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	int j,last,first,reject,test,ta,tb;
	real_t xt,h,hmin,ind,hl,ht,absh,fhm,discr,tol,mu,mu1,fh,e1,e2,
			*k0,*k1,*k2,*k3,*k4;

	k0=allocate_real_vector(1,n);
	k1=allocate_real_vector(1,n);
	k2=allocate_real_vector(1,n);
	k3=allocate_real_vector(1,n);
	k4=allocate_real_vector(1,n);
	if (fi) {
		data[3]=(*xe)-(*x);
		data[4]=data[5]=data[6]=0.0;
	}
	absh=h=fabs(data[3]);
	if (*xe < *x) h = -h;
	ind=fabs((*xe)-(*x));
	hmin=ind*data[1]+data[2];
	e1=12.0*data[1]/ind;
	e2=12.0*data[2]/ind;
	first=1;
	reject=0;
	test=1;
	if (fi) {
		last=1;
		test=0;
	}
	while (1) {
		if (test) {
			absh=fabs(h);
			if (absh < hmin) {
				h = (*xe == *x) ? 0.0 : ((*xe > *x) ? hmin : -hmin);
				absh=hmin;
			}
			ta=(h >= (*xe)-(*x));
			tb=(h >= 0.0);
			if ((ta && tb) || (!(ta || tb))) {
				last=1;
				h=(*xe)-(*x);
				absh=fabs(h);
			} else
				last=0;
		}
		test=1;
		if (!reject) {
			for (j=1; j<=n; j++) k0[j]=y[j];
			(*der)(n,*x,k0);
		}
		ht=0.184262134833347*h;
		xt = *x+ht;
		for (j=1; j<=n; j++) k1[j]=k0[j]*ht+y[j];
		(*der)(n,xt,k1);
		ht=0.690983005625053e-1*h;
		xt=4.0*ht+(*x);
		for (j=1; j<=n; j++) k2[j]=(3.0*k1[j]+k0[j])*ht+y[j];
		(*der)(n,xt,k2);
		xt=0.5*h+(*x);
		ht=0.1875*h;
		for (j=1; j<=n; j++)
			k3[j]=((1.74535599249993*k2[j]-k1[j])*2.23606797749979+
						k0[j])*ht+y[j];
		(*der)(n,xt,k3);
		xt=0.723606797749979*h+(*x);
		ht=0.4*h;
		for (j=1; j<=n; j++)
			k4[j]=(((0.517595468166681*k0[j]-k1[j])*0.927050983124840+
						k2[j])*1.46352549156242+k3[j])*ht+y[j];
		(*der)(n,xt,k4);
		xt = (last ? *xe : *x+h);
		ht=2.0*h;
		for (j=1; j<=n; j++)
			k1[j]=((((2.0*k4[j]+k2[j])*0.412022659166595+k1[j])*
						2.23606797749979-k0[j])*0.375-k3[j])*ht+y[j];
		(*der)(n,xt,k1);
		reject=0;
		fhm=0.0;
		for (j=1; j<=n; j++) {
			discr=fabs((1.6*k3[j]-k2[j]-k4[j])*5.0+k0[j]+k1[j]);
			tol=fabs(k0[j])*e1+e2;
			reject = (discr > tol || reject);
			fh=discr/tol;
			if (fh > fhm) fhm=fh;
		}
		mu=1.0/(1.0+fhm)+0.45;
		if (reject) {
			data[5] += 1.0;
			if (absh <= hmin) {
				data[6] += 1.0;
				hl=h;
				reject=0;
				first=1;
				data[3]=hl;
				data[4] += 1.0;
				*x=xt;
				(*out)(n,*x,*xe,y,data);
				if (*x == *xe) break;
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
			ht=hl/12.0;
			for (j=1; j<=n; j++)
				y[j]=((k2[j]+k4[j])*5.0+k0[j]+k1[j])*ht+y[j];
			data[3]=hl;
			data[4] += 1.0;
			*x=xt;
			(*out)(n,*x,*xe,y,data);
			if (*x == *xe) break;
		}
	}
	free_real_vector(k0,1);
	free_real_vector(k1,1);
	free_real_vector(k2,1);
	free_real_vector(k3,1);
	free_real_vector(k4,1);
}
