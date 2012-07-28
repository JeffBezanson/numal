#include "../real.h"


void reaveches(real_t **a, int n, real_t lambda, real_t em[], real_t v[])
{
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	int i,i1,j,count,max,*p;
	real_t m,r,norm,machtol,tol;

	p=allocate_integer_vector(1,n);
	norm=em[1];
	machtol=em[0]*norm;
	tol=em[6]*norm;
	max=em[8];
	a[1][1] -= lambda;
	for (i=1; i<=n-1; i++) {
		i1=i+1;
		r=a[i][i];
		m=a[i1][i];
		if (fabs(m) < machtol) m=machtol;
		p[i] = (fabs(m) <= fabs(r));
		if (p[i]) {
			a[i1][i] = m /= r;
			for (j=i1; j<=n; j++)
				a[i1][j]=((j > i1) ? a[i1][j] : a[i1][j]-lambda)-m*a[i][j];
		} else {
			a[i][i]=m;
			a[i1][i] = m = r/m;
			for (j=i1; j<=n; j++) {
				r = (j > i1) ? a[i1][j] : a[i1][j]-lambda;
				a[i1][j]=a[i][j]-m*r;
				a[i][j]=r;
			}
		}
	}
	if (fabs(a[n][n]) < machtol) a[n][n]=machtol;
	for (j=1; j<=n; j++) v[j]=1.0;
	count=0;
	do {
		count++;
		if (count > max) break;
		for (i=1; i<=n-1; i++) {
			i1=i+1;
			if (p[i])
				v[i1] -= a[i1][i]*v[i];
			else {
				r=v[i1];
				v[i1]=v[i]-a[i1][i]*r;
				v[i]=r;
			}
		}
		for (i=n; i>=1; i--)
			v[i]=(v[i]-matvec(i+1,n,i,a,v))/a[i][i];
		r=1.0/sqrt(vecvec(1,n,0,v,v));
		for (j=1; j<=n; j++) v[j] *= r;
	} while (r > tol);
	em[7]=r;
	em[9]=count;
	free_integer_vector(p,1);
}

