#include "../real.h"


void reccof(int n, int m, real_t *x, real_t (*wx)(real_t), real_t b[],
				real_t c[], real_t l[], int sym)
{
	real_t ortpol(int, real_t, real_t [], real_t []);
	int i,j,up;
	real_t r,s,pim,h,hh,arg,sa,temp;

	pim=4.0*atan(1.0)/m;
	if (sym) {
		for (j=0; j<=n; j++) {
			r=b[j]=0.0;
			up=m/2;
			for (i=1; i<=up; i++) {
				arg=(i-0.5)*pim;
				*x=cos(arg);
				temp=ortpol(j,*x,b,c);
				r += sin(arg)*(*wx)(*x)*temp*temp;
			}
			if (up*2 == m)
				l[j]=2.0*r*pim;
			else {
				*x=0.0;
				temp=ortpol(j,0.0,b,c);
				l[j]=(2.0*r+(*wx)(*x)*temp*temp)*pim;
			}
			c[j] = (j == 0) ? 0.0 : l[j]/l[j-1];
		}
	} else
		for (j=0; j<=n; j++) {
			r=s=0.0;
			up=m/2;
			for (i=1; i<=up; i++) {
				arg=(i-0.5)*pim;
				sa=sin(arg);
				*x=cos(arg);
				temp=ortpol(j,*x,b,c);
				h=(*wx)(*x)*temp*temp;
				*x = -(*x);
				temp=ortpol(j,*x,b,c);
				hh=(*wx)(*x)*temp*temp;
				r += (h+hh)*sa;
				s += (hh-h)*(*x)*sa;
			}
			b[j]=s*pim;
			if (up*2 == m)
				l[j]=r*pim;
			else {
				*x=0.0;
				temp=ortpol(j,0.0,b,c);
				l[j]=(r+(*wx)(*x)*temp*temp)*pim;
			}
			c[j] = (j == 0) ? 0.0 : l[j]/l[j-1];
		}
}
