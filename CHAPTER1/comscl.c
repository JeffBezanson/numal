#include "../real.h"


void comscl(real_t **a, int n, int n1, int n2, real_t im[])
{
	int i,j,k;
	real_t s,u,v,w,aij,aij1;

	for (j=n1; j<=n2; j++) {
		s=0.0;
		if (im[j] != 0.0) {
			for (i=1; i<=n; i++) {
				aij=a[i][j];
				aij1=a[i][j+1];
				u=aij*aij+aij1*aij1;
				if (u > s) {
					s=u;
					k=i;
				}
			}
			if (s != 0.0) {
				v=a[k][j]/s;
				w = -a[k][j+1]/s;
				for (i=1; i<=n; i++) {
					u=a[i][j];
					s=a[i][j+1];
					a[i][j]=u*v-s*w;
					a[i][j+1]=u*w+s*v;
				}
			}
			j++;
		} else {
			for (i=1; i<=n; i++)
				if (fabs(a[i][j]) > fabs(s)) s=a[i][j];
			if (s != 0.0)
				for (i=1; i<=n; i++)
					a[i][j] /= s;
		}
	}
}
