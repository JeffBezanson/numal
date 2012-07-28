#include "../real.h"


int reavalqri(real_t **a, int n, real_t em[], real_t val[])
{
	void rotcol(int, int, int, int, real_t **, real_t, real_t);
	void rotrow(int, int, int, int, real_t **, real_t, real_t);
	int n1,i,i1,q,max,count;
	real_t det,w,shift,kappa,nu,mu,r,tol,delta,machtol,s;

	machtol=em[0]*em[1];
	tol=em[1]*em[2];
	max=em[4];
	count=0;
	r=0.0;
	do {
		n1=n-1;
		i=n;
		do{
			q=i;
			i--;
		} while ((i >= 1) ? (fabs(a[i+1][i]) > tol) : 0);
		if (q > 1)
			if (fabs(a[q][q-1]) > r) r=fabs(a[q][q-1]);
		if (q == n) {
			val[n]=a[n][n];
			n=n1;
		} else {
			delta=a[n][n]-a[n1][n1];
			det=a[n][n1]*a[n1][n];
			if (fabs(delta) < machtol)
				s=sqrt(det);
			else {
				w=2.0/delta;
				s=w*w*det+1.0;
				s = (s <= 0.0) ? -delta*0.5 : w*det/(sqrt(s)+1.0);
			}
			if (q == n1) {
				val[n]=a[n][n]+s;
				val[n1]=a[n1][n1]-s;
				n -= 2;
			} else {
				count++;
				if (count > max) break;
				shift=a[n][n]+s;
				if (fabs(delta) < tol) {
					w=a[n1][n1]-s;
					if (fabs(w) < fabs(shift)) shift=w;
				}
				a[q][q] -= shift;
				for (i=q; i<=n-1; i++) {
					i1=i+1;
					a[i1][i1] -= shift;
					kappa=sqrt(a[i][i]*a[i][i]+a[i1][i]*a[i1][i]);
					if (i > q) {
						a[i][i-1]=kappa*nu;
						w=kappa*mu;
					} else
						w=kappa;
					mu=a[i][i]/kappa;
					nu=a[i1][i]/kappa;
					a[i][i]=w;
					rotrow(i1,n,i,i1,a,mu,nu);
					rotcol(q,i,i,i1,a,mu,nu);
					a[i][i] += shift;
				}
				a[n][n-1]=a[n][n]*nu;
				a[n][n]=a[n][n]*mu+shift;
			}
		}
	} while (n > 0);
	em[3]=r;
	em[5]=count;
	return n;
}

