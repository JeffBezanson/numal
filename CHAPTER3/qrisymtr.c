#include "../real.h"


int qrisymtri(real_t **a, int n, real_t d[], real_t b[], real_t bb[],
					real_t em[])
{
	void rotcol(int, int, int, int, real_t **, real_t, real_t);
	int j,j1,k,m,m1,count,max;
	real_t bbmax,r,s,sin,t,cos,oldcos,g,p,w,tol,tol2,lambda,dk1;

	tol=em[2]*em[1];
	tol2=tol*tol;
	count=0;
	bbmax=0.0;
	max=em[4];
	m=n;
	do {
		k=m;
		m1=m-1;
		while (1) {
			k--;
			if (k <= 0) break;
			if (bb[k] < tol2) {
				if (bb[k] > bbmax) bbmax=bb[k];
				break;
			}
		}
		if (k == m1)
			m=m1;
		else {
			t=d[m]-d[m1];
			r=bb[m1];
			if (fabs(t) < tol)
				s=sqrt(r);
			else {
				w=2.0/t;
				s=w*r/(sqrt(w*w*r+1.0)+1.0);
			}
			if (k == m-2) {
				d[m] += s;
				d[m1] -= s;
				t = -s/b[m1];
				r=sqrt(t*t+1.0);
				cos=1.0/r;
				sin=t/r;
				rotcol(1,n,m1,m,a,cos,sin);
				m -= 2;
			} else {
				count++;
				if (count > max) break;
				lambda=d[m]+s;
				if (fabs(t) < tol) {
					w=d[m1]-s;
					if (fabs(w) < fabs(lambda)) lambda=w;
				}
				k++;
				t=d[k]-lambda;
				cos=1.0;
				w=b[k];
				p=sqrt(t*t+w*w);
				j1=k;
				for (j=k+1; j<=m; j++) {
					oldcos=cos;
					cos=t/p;
					sin=w/p;
					dk1=d[j]-lambda;
					t *= oldcos;
					d[j1]=(t+dk1)*sin*sin+lambda+t;
					t=cos*dk1-sin*w*oldcos;
					w=b[j];
					p=sqrt(t*t+w*w);
					g=b[j1]=sin*p;
					bb[j1]=g*g;
					rotcol(1,n,j1,j,a,cos,sin);
					j1=j;
				}
				d[m]=cos*t+lambda;
				if (t < 0.0) b[m1] = -g;
			}
		}
	} while (m > 0);
	em[3]=sqrt(bbmax);
	em[5]=count;
	return m;
}

