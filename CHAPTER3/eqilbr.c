#include "../real.h"


void eqilbr(real_t **a, int n, real_t em[], real_t d[], int inter[])
{
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t mattam(int, int, int, int, real_t **, real_t **);
	void ichcol(int, int, int, int, real_t **);
	void ichrow(int, int, int, int, real_t **);
	int i,im,i1,p,q,j,t,count,exponent,ni;
	real_t c,r,eps,omega,factor,di;

	factor=1.0/(2.0*log(2.0));
	eps=em[0];
	omega=1.0/eps;
	t=p=1;
	q=ni=i=n;
	count=((n+1)*n)/2;
	for (j=1; j<=n; j++) {
		d[j]=1.0;
		inter[j]=0;
	}
	i = (i < q) ? i+1 : p;
	while (count > 0 && ni > 0) {
		count--;
		im=i-1;
		i1=i+1;
		c=sqrt(tammat(p,im,i,i,a,a)+tammat(i1,q,i,i,a,a));
		r=sqrt(mattam(p,im,i,i,a,a)+mattam(i1,q,i,i,a,a));
		if (c*omega <= r*eps) {
			inter[t]=i;
			ni=q-p;
			t++;
			if (p != i) {
				ichcol(1,n,p,i,a);
				ichrow(1,n,p,i,a);
				di=d[i];
				d[i]=d[p];
				d[p]=di;
			}
			p++;
		} else
			if (r*omega <= c*eps) {
				inter[t] = -i;
				ni=q-p;
				t++;
				if (q != i) {
					ichcol(1,n,q,i,a);
					ichrow(1,n,q,i,a);
					di=d[i];
					d[i]=d[q];
					d[q]=di;
				}
				q--;
			} else {
				exponent=log(r/c)*factor;
				if (fabs(exponent) > 1.0) {
					ni=q-p;
					c=pow(2.0,exponent);
					r=1.0/c;
					d[i] *= c;
					for (j=1; j<=im; j++) {
						a[j][i] *=c;
						a[i][j] *= r;
					}
					for (j=i1; j<=n; j++) {
						a[j][i] *=c;
						a[i][j] *= r;
					}
				} else
					ni--;
			}
		i = (i < q) ? i+1 : p;
	}
}
