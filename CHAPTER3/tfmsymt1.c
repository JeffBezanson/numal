#include "../real.h"


void tfmsymtri1(real_t a[], int n, real_t d[], real_t b[], real_t bb[],
					real_t em[])
{
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t seqvec(int, int, int, int, real_t [], real_t []);
	void elmvec(int, int, int, real_t [], real_t [], real_t);
	int i,j,r,r1,p,q,ti,tj;
	real_t s,w,x,a1,b0,bb0,norm,machtol;

	norm=0.0;
	tj=0;
	for (j=1; j<=n; j++) {
		w=0.0;
		for (i=1; i<=j; i++) w += fabs(a[i+tj]);
		tj += j;
		ti=tj+j;
		for (i=j+1; i<=n; i++) {
			w += fabs(a[ti]);
			ti += i;
		}
		if (w > norm) norm=w;
	}
	machtol=em[0]*norm;
	em[1]=norm;
	q=((n+1)*n)/2;
	r=n;
	for (r1=n-1; r1>=1; r1--) {
		p=q-r;
		d[r]=a[q];
		x=vecvec(p+1,q-2,0,a,a);
		a1=a[q-1];
		if (sqrt(x) <= machtol) {
			b0=b[r1]=a1;
			bb[r1]=b0*b0;
			a[q]=1.0;
		} else {
			bb0=bb[r1]=a1*a1+x;
			b0 = (a1 > 0.0) ? -sqrt(bb0) : sqrt(bb0);
			a1=a[q-1]=a1-b0;
			w=a[q]=1.0/(a1*b0);
			tj=0;
			for (j=1; j<=r1; j++) {
				ti=tj+j;
				s=vecvec(tj+1,ti,p-tj,a,a);
				tj=ti+j;
				b[j]=(seqvec(j+1,r1,tj,p,a,a)+s)*w;
				tj=ti;
			}
			elmvec(1,r1,p,b,a,vecvec(1,r1,p,b,a)*w*0.5);
			tj=0;
			for (j=1; j<=r1; j++) {
				ti=tj+j;
				elmvec(tj+1,ti,p-tj,a,a,b[j]);
				elmvec(tj+1,ti,-tj,a,b,a[j+p]);
				tj=ti;
			}
			b[r1]=b0;
		}
		q=p;
		r=r1;
	}
	d[1]=a[1];
	a[1]=1.0;
	b[n]=bb[n]=0.0;
}
