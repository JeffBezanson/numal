#include "../real.h"


int qrivalsymtri(real_t d[], real_t bb[], int n, real_t em[])
{
	int i,i1,low,oldlow,n1,count,max;
	real_t bbtol,bbmax,bbi,bbn1,machtol,dn,delta,f,num,shift,g,h,
			t,p,r,s,c,oldg;

	t=em[2]*em[1];
	bbtol=t*t;
	machtol=em[0]*em[1];
	max=em[4];
	bbmax=0.0;
	count=0;
	oldlow=n;
	n1=n-1;
	while (n > 0) {
		i=n;
		do {
			low=i;
			i--;
		} while ((i >= 1) ? bb[i] > bbtol : 0);
		if (low > 1)
			if (bb[low-1] > bbmax) bbmax=bb[low-1];
		if (low == n)
			n=n1;
		else {
			dn=d[n];
			delta=d[n1]-dn;
			bbn1=bb[n1];
			if (fabs(delta) < machtol)
				r=sqrt(bbn1);
			else {
				f=2.0/delta;
				num=bbn1*f;
				r = -num/(sqrt(num*f+1.0)+1.0);
			}
			if (low == n1) {
				d[n]=dn+r;
				d[n1] -= r;
				n -= 2;
			} else {
				count++;
				if (count > max) break;
				if (low < oldlow) {
					shift=0.0;
					oldlow=low;
				} else
					shift=dn+r;
				h=d[low]-shift;
				if (fabs(h) < machtol)
					h = (h <= 0.0) ? -machtol : machtol;
				g=h;
				t=g*h;
				bbi=bb[low];
				p=t+bbi;
				i1=low;
				for (i=low+1; i<=n; i++) {
					s=bbi/p;
					c=t/p;
					h=d[i]-shift-bbi/h;
					if (fabs(h) < machtol)
						h = (h <= 0.0) ? -machtol : machtol;
					oldg=g;
					g=h*c;
					t=g*h;
					d[i1]=oldg-g+d[i];
					bbi = (i == n) ? 0.0 : bb[i];
					p=t+bbi;
					bb[i1]=s*p;
					i1=i;
				}
				d[n]=g+shift;
			}
		}
		n1=n-1;
	}
	em[3]=sqrt(bbmax);
	em[5]=count;
	return n;
}

