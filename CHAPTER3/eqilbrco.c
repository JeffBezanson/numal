#include "../real.h"
#include <stdlib.h>


void eqilbrcom(real_t **a1, real_t **a2, int n, real_t em[], real_t d[],
					int inter[])
{
	void ichcol(int, int, int, int, real_t **);
	void ichrow(int, int, int, int, real_t **);
	real_t tammat(int, int, int, int, real_t **, real_t **);
	real_t mattam(int, int, int, int, real_t **, real_t **);
	int i,p,q,j,t,count,exponent,ni,im,i1;
	real_t c,r,eps,di;

	eps=em[0]*em[0];
	t=p=1;
	q=ni=i=n;
	count=em[6];
	for (j=1; j<=n; j++) {
		d[j]=1.0;
		inter[j]=0;
	}
	i = (i < q) ? i+1 : p;
	while (count > 0 && ni > 0) {
		count--;
		im=i-1;
		i1=i+1;
		c=tammat(p,im,i,i,a1,a1)+tammat(i1,q,i,i,a1,a1)+
			tammat(p,im,i,i,a2,a2)+tammat(i1,q,i,i,a2,a2);
		r=mattam(p,im,i,i,a1,a1)+mattam(i1,q,i,i,a1,a1)+
			mattam(p,im,i,i,a2,a2)+mattam(i1,q,i,i,a2,a2);
		if (c/eps <= r) {
			inter[t]=i;
			ni=q-p;
			t++;
			if (p != i) {
				ichcol(1,n,p,i,a1);
				ichrow(1,n,p,i,a1);
				ichcol(1,n,p,i,a2);
				ichrow(1,n,p,i,a2);
				di=d[i];
				d[i]=d[p];
				d[p]=di;
			}
			p++;
		} else
			if (r/eps <= c) {
				inter[t] = -i;
				ni=q-p;
				t++;
				if (q != i) {
					ichcol(1,n,q,i,a1);
					ichrow(1,n,q,i,a1);
					ichcol(1,n,q,i,a2);
					ichrow(1,n,q,i,a2);
					di=d[i];
					d[i]=d[q];
					d[q]=di;
				}
				q--;
			} else {
				exponent=ceil(log(r/c)*0.36067);
				if (abs(exponent) > 1) {
					ni=q-p;
					c=pow(2.0,exponent);
					d[i] *= c;
					for (j=1; j<=im; j++) {
						a1[j][i] *= c;
						a1[i][j] /= c;
						a2[j][i] *= c;
						a2[i][j] /= c;
					}
					for (j=i1; j<=n; j++) {
						a1[j][i] *= c;
						a1[i][j] /= c;
						a2[j][i] *= c;
						a2[i][j] /= c;
					}
				} else
					ni--;
			}
		i = (i < q) ? i+1 : p;
	}
	em[7]=em[6]-count;
}
