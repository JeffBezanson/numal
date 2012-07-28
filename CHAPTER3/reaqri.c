#include "../real.h"


int reaqri(real_t **a, int n, real_t em[], real_t val[], real_t **vec)
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	real_t matvec(int, int, int, real_t **, real_t []);
	void rotcol(int, int, int, int, real_t **, real_t, real_t);
	void rotrow(int, int, int, int, real_t **, real_t, real_t);
	int m1,i,i1,m,j,q,max,count;
	real_t w,shift,kappa,nu,mu,r,tol,s,machtol,elmax,t,delta,det,*tf;

	tf=allocate_real_vector(1,n);
	machtol=em[0]*em[1];
	tol=em[1]*em[2];
	max=em[4];
	count=0;
	elmax=0.0;
	m=n;
	for (i=1; i<=n; i++) {
		vec[i][i]=1.0;
		for (j=i+1; j<=n; j++) vec[i][j]=vec[j][i]=0.0;
	}
	do {
		m1=m-1;
		i=m;
		do {
			q=i;
			i--;
		} while ((i >= 1) ? (fabs(a[i+1][i]) > tol) : 0);
		if (q > 1)
			if (fabs(a[q][q-1]) > elmax) elmax=fabs(a[q][q-1]);
		if (q == m) {
			val[m]=a[m][m];
			m=m1;
		} else {
			delta=a[m][m]-a[m1][m1];
			det=a[m][m1]*a[m1][m];
			if (fabs(delta) < machtol)
				s=sqrt(det);
			else {
				w=2.0/delta;
				s=w*w*det+1.0;
				s = (s <= 0.0) ? -delta*0.5 : w*det/(sqrt(s)+1.0);
			}
			if (q == m1) {
				val[m] = a[m][m] += s;
				val[q] = a[q][q] -= s;
				t = (fabs(s) < machtol) ? (s+delta)/a[m][q] : a[q][m]/s;
				r=sqrt(t*t+1.0);
				nu=1.0/r;
				mu = -t*nu;
				a[q][m] -= a[m][q];
				rotrow(q+2,n,q,m,a,mu,nu);
				rotcol(1,q-1,q,m,a,mu,nu);
				rotcol(1,n,q,m,vec,mu,nu);
				m -= 2;
			} else {
				count++;
				if (count > max) {
					em[3]=elmax;
					em[5]=count;
					free_real_vector(tf,1);
					return m;
				}
				shift=a[m][m]+s;
				if (fabs(delta) < tol) {
					w=a[m1][m1]-s;
					if (fabs(w) < fabs(shift)) shift=w;
				}
				a[q][q] -= shift;
				for (i=q; i<=m1; i++) {
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
					rotcol(1,i,i,i1,a,mu,nu);
					a[i][i] += shift;
					rotcol(1,n,i,i1,vec,mu,nu);
				}
				a[m][m1]=a[m][m]*nu;
				a[m][m]=a[m][m]*mu+shift;
			}
		}
	} while (m > 0);
	for (j=n; j>=2; j--) {
		tf[j]=1.0;
		t=a[j][j];
		for (i=j-1; i>=1; i--) {
			delta=t-a[i][i];
			tf[i]=matvec(i+1,j,i,a,tf)/
						((fabs(delta) < machtol) ? machtol : delta);
		}
		for (i=1; i<=n; i++) vec[i][j]=matvec(1,j,i,vec,tf);
	}
	em[3]=elmax;
	em[5]=count;
	free_real_vector(tf,1);
	return m;
}

