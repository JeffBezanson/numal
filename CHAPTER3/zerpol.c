#include "../real.h"


real_t *allocate_real_vector(int, int);
void free_real_vector(real_t *, int);
real_t comabs(real_t, real_t);
void comsqrt(real_t, real_t, real_t *, real_t *);
int zerpolfunction(int, real_t [], real_t [], real_t,
				   real_t, real_t, int *, real_t *);

int zerpol(int n, real_t a[], real_t em[], real_t re[], real_t im[],
				real_t d[])
{
	int i,totit,it,fail,start,up,max,giex,itmax,control,ih,m,split;
	real_t x,y,newf,oldf,maxrad,ae,tol,h1,h2,ln2,f[6],tries[11],
			h,side,s1re,s1im,s2re,s2im,dx,dy,h3,h4,h5,h6,*b;

	totit=it=fail=up=start=0;
	ln2=log(2.0);
	newf=FLT_MAX;
	ae=FLT_MIN;
	giex=log(newf)/ln2-40.0;
	tol=em[0];
	itmax=em[1];
	for (i=0; i<=n; i++) d[i]=a[n-i];
	if (n <= 0)
		fail=1;
	else
		if (d[0] == 0.0) fail=2;
	if (fail > 0) {
		em[2]=fail;
		em[3]=start;
		em[4]=totit;
		for (i=(n-1)/2; i>=0; i--) {
			tol=d[i];
			d[i]=d[n-i];
			d[n-i]=tol;
		}
		return n;
	}
	while (d[n] == 0.0 && n > 0) {
		re[n]=im[n]=0.0;
		n--;
	}
	x=y=0.0;
	while (n > 2) {
		/* control */
		if (it > itmax) {
			totit += it;
			fail=3;
			em[2]=fail;
			em[3]=start;
			em[4]=totit;
			for (i=(n-1)/2; i>=0; i--) {
				tol=d[i];
				d[i]=d[n-i];
				d[n-i]=tol;
			}
			return n;
		} else
			if (it == 0) {
				maxrad=0.0;
				max=(giex-log(fabs(d[0]))/ln2)/n;
				for (i=1; i<=n; i++) {
					h1 = (d[i] == 0.0) ? 0.0 : exp(log(fabs(d[i]/d[0]))/i);
					if (h1 > maxrad) maxrad=h1;
				}
				for (i=1; i<=n-1; i++)
					if (d[i] != 0.0) {
						ih=(giex-log(fabs(d[i]))/ln2)/(n-i);
						if (ih < max) max=ih;
					}
				max=max*ln2/log(n);
				side = -d[1]/d[0];
				side = (fabs(side) < tol) ? 0.0 :
							((side > 0.0) ? 1.0 : -1.0);
				if (side == 0.0) {
					tries[7]=tries[2]=maxrad;
					tries[9] = -maxrad;
					tries[6]=tries[4]=tries[3]=maxrad/sqrt(2.0);
					tries[5] = -tries[3];
					tries[10]=tries[8]=tries[1]=0.0;
				} else {
					tries[8]=tries[4]=maxrad/sqrt(2.0);
					tries[1]=side*maxrad;
					tries[3]=tries[4]*side;
					tries[6]=maxrad;
					tries[7] = -tries[3];
					tries[9] = -tries[1];
					tries[2]=tries[5]=tries[10]=0.0;
				}
				if (comabs(x,y) > 2.0*maxrad) x=y=0.0;
				control=0;
			} else {
				if (it > 1 && newf >= oldf) {
					up++;
					if (up == 5 && start < 5) {
						start++;
						up=0;
						x=tries[2*start-1];
						y=tries[2*start];
						control=0;
					} else
						control=1;
				} else
					control=1;
			}	/* end of control */
		if (control) {
			/* laguerre */
			if (fabs(f[0]) > fabs(f[1])) {
				h1=f[0];
				h6=f[1]/h1;
				h2=f[2]+h6*f[3];
				h3=f[3]-h6*f[2];
				h4=f[4]+h6*f[5];
				h5=f[5]-h6*f[4];
				h6=h6*f[1]+h1;
			} else {
				h1=f[1];
				h6=f[0]/h1;
				h2=h6*f[2]+f[3];
				h3=h6*f[3]-f[2];
				h4=h6*f[4]+f[5];
				h5=h6*f[5]-f[4];
				h6=h6*f[0]+f[1];
			}
			s1re=h2/h6;
			s1im=h3/h6;
			h2=s1re*s1re-s1im*s1im;
			h3=2.0*s1re*s1im;
			s2re=h2-h4/h6;
			s2im=h3-h5/h6;
			h1=s2re*s2re+s2im*s2im;
			h1 = (h1 != 0.0) ? (s2re*h2+s2im*h3)/h1 : 1.0;
			m = (h1 > n-1) ? ((n > 1) ? n-1 : 1) : ((h1 > 1.0) ? h1 : 1);
			h1=(real_t)(n-m)/(real_t) m;
			comsqrt(h1*(n*s2re-h2),h1*(n*s2im-h3),&h2,&h3);
			if (s1re*h2+s1im*h3 < 0.0) {
				h2 = -h2;
				h3 = -h3;
			}
			h2 += s1re;
			h3 += s1im;
			h1=h2*h2+h3*h3;
			if (h1 == 0.0) {
				dx = -n;
				dy=n;
			} else {
				dx = -n*h2/h1;
				dy=n*h3/h1;
			}
			h1=fabs(x)*tol+ae;
			h2=fabs(y)*tol+ae;
			if (fabs(dx) < h1 && fabs(dy) < h2) {
				dx = (dx == 0.0) ? h1 : ((dx > 0.0) ? h1 : -h1);
				dy = (dy == 0.0) ? h2 : ((dy > 0.0) ? h2 : -h2);
			}
			x += dx;
			y += dy;
			if (comabs(x,y) > 2.0*maxrad) {
				h1 = (fabs(x) > fabs(y)) ? fabs(x) : fabs(y);
				h2=log(h1)/ln2+1.0-max;
				if (h2 > 0.0) {
					h2=pow(2.0,h2);
					x /= h2;
					y /= h2;
				}
			}	/* end of laguerre */
		}
		oldf=newf;
		if (zerpolfunction(n,d,f,x,y,tol,&it,&newf)) {
			if (y != 0.0 && fabs(y) < 0.1) {
				h=y;
				y=0.0;
				if (!zerpolfunction(n,d,f,x,y,tol,&it,&newf)) y=h;
			}
			re[n]=x;
			im[n]=y;
			if (y != 0.0) {
				re[n-1]=x;
				im[n-1] = -y;
			}
			/* deflation */
			if (x == 0.0 && y == 0.0)
				n--;
			else {
				b=allocate_real_vector(0,n-1);
				if (y == 0.0) {
					n--;
					b[n] = -d[n+1]/x;
					for (i=1; i<=n; i++) b[n-i]=(b[n-i+1]-d[n-i+1])/x;
					for (i=1; i<=n; i++) d[i] += d[i-1]*x;
				} else {
					h1 = -2.0*x;
					h2=x*x+y*y;
					n -= 2;
					b[n]=d[n+2]/h2;
					b[n-1]=(d[n+1]-h1*b[n])/h2;
					for (i=2; i<=n; i++)
						b[n-i]=(d[n-i+2]-h1*b[n-i+1]-b[n-i+2])/h2;
					d[1] -= h1*d[0];
					for (i=2; i<=n; i++) d[i] -= h1*d[i-1]+h2*d[i-2];
				}
				split=n;
				h2=fabs(d[n]-b[n])/(fabs(d[n])+fabs(b[n]));
				for (i=n-1; i>=0; i--) {
					h1=fabs(d[i])+fabs(b[i]);
					if (h1 > tol) {
						h1=fabs(d[i]-b[i])/h1;
						if (h1 < h2) {
							h2=h1;
							split=i;
						}
					}
				}
				for (i=split+1; i<=n; i++) d[i]=b[i];
				d[split]=(d[split]+b[split])/2.0;
				free_real_vector(b,0);
			}	/* end of deflation */
			totit += it;
			up=start=it=0;
		}
	}
	if (n == 1) {
		re[1] = -d[1]/d[0];
		im[1]=0.0;
	} else {
		h1 = -0.5*d[1]/d[0];
		h2=h1*h1-d[2]/d[0];
		if (h2 >= 0.0) {
			re[2] = (h1 < 0.0) ? h1-sqrt(h2) : h1+sqrt(h2);
			re[1]=d[2]/(d[0]*re[2]);
			im[2]=im[1]=0.0;
		} else {
			re[2]=re[1]=h1;
			im[2]=sqrt(-h2);
			im[1] = -im[2];
		}
	}
	em[2]=fail;
	em[3]=start;
	em[4]=totit;
	return 0;
}

int zerpolfunction(int n, real_t d[], real_t f[], real_t x,
						real_t y, real_t tol, int *it, real_t *newf)
{
	/* this function is used internally by ZERPOL */

	int k,m1,m2;
	real_t p,q,qsqrt,f01,f02,f03,f11,f12,f13,f21,f22,f23,stop;

	(*it)++;
	p=2.0*x;
	q = -(x*x+y*y);
	qsqrt=sqrt(-q);
	f01=f11=f21=d[0];
	f02=f12=f22=0.0;
	m1=n-4;
	m2=n-2;
	stop=fabs(f01)*0.8;
	for (k=1; k<=m1; k++) {
		f03=f02;
		f02=f01;
		f01=d[k]+p*f02+q*f03;
		f13=f12;
		f12=f11;
		f11=f01+p*f12+q*f13;
		f23=f22;
		f22=f21;
		f21=f11+p*f22+q*f23;
		stop=qsqrt*stop+fabs(f01);
	}
	if (m1 < 0) m1=0;
	for (k=m1+1; k<=m2; k++) {
		f03=f02;
		f02=f01;
		f01=d[k]+p*f02+q*f03;
		f13=f12;
		f12=f11;
		f11=f01+p*f12+q*f13;
		stop=qsqrt*stop+fabs(f01);
	}
	if (n == 3) f21=0.0;
	f03=f02;
	f02=f01;
	f01=d[n-1]+p*f02+q*f03;
	f[0]=d[n]+x*f01+q*f02;
	f[1]=y*f01;
	f[2]=f01-2.0*f12*y*y;
	f[3]=2.0*y*(-x*f12+f11);
	f[4]=2.0*(-x*f12+f11)-8.0*y*y*(-x*f22+f21);
	f[5]=y*(6.0*f12-8.0*y*y*f22);
	stop=qsqrt*(qsqrt*stop+fabs(f01))+fabs(f[0]);
	*newf=f02=comabs(f[0],f[1]);
	return (f02 < (2.0*fabs(x*f01)-8.0*(fabs(f[0])+fabs(f01)*qsqrt)+
				10.0*stop)*tol*pow(1.0+tol,4*n+3.0));
}

