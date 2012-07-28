#include "../real.h"


void comveches(real_t **a, int n, real_t lambda, real_t mu,
					real_t em[], real_t u[], real_t v[])
{
	real_t *allocate_real_vector(int, int);
	int *allocate_integer_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	real_t matvec(int, int, int, real_t **, real_t []);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int i,i1,j,count,max,*p;
	real_t aa,bb,d,m,r,s,w,x,y,norm,machtol,tol,*g,*f;

	p=allocate_integer_vector(1,n);
	g=allocate_real_vector(1,n);
	f=allocate_real_vector(1,n);
	norm=em[1];
	machtol=em[0]*norm;
	tol=em[6]*norm;
	max=em[8];
	for (i=2; i<=n; i++) {
		f[i-1]=a[i][i-1];
		a[i][1]=0.0;
	}
	aa=a[1][1]-lambda;
	bb = -mu;
	for (i=1; i<=n-1; i++) {
		i1=i+1;
		m=f[i];
		if (fabs(m) < machtol) m=machtol;
		a[i][i]=m;
		d=aa*aa+bb*bb;
		p[i] = (fabs(m) < sqrt(d));
		if (p[i]) {
			f[i]=r=m*aa/d;
			g[i] = s = -m*bb/d;
			w=a[i1][i];
			x=a[i][i1];
			a[i1][i]=y=x*s+w*r;
			a[i][i1]=x=x*r-w*s;
			aa=a[i1][i1]-lambda-x;
			bb = -(mu+y);
			for (j=i+2; j<=n; j++) {
				w=a[j][i];
				x=a[i][j];
				a[j][i]=y=x*s+w*r;
				a[i][j]=x=x*r-w*s;
				a[j][i1] = -y;
				a[i1][j] -= x;
			}
		} else {
			f[i]=r=aa/m;
			g[i]=s=bb/m;
			w=a[i1][i1]-lambda;
			aa=a[i][i1]-r*w-s*mu;
			a[i][i1]=w;
			bb=a[i1][i]-s*w+r*mu;
			a[i1][i] = -mu;
			for (j=i+2; j<=n; j++) {
				w=a[i1][j];
				a[i1][j]=a[i][j]-r*w;
				a[i][j]=w;
				a[j][i1]=a[j][i]-s*w;
				a[j][i]=0.0;
			}
		}
	}
	p[n]=1;
	d=aa*aa+bb*bb;
	if (d < machtol*machtol) {
		aa=machtol;
		bb=0.0;
		d=machtol*machtol;
	}
	a[n][n]=d;
	f[n]=aa;
	g[n] = -bb;
	for (i=1; i<=n; i++) {
		u[i]=1.0;
		v[i]=0.0;
	}
	count=0;
	do {
		if (count > max) break;
		for (i=1; i<=n; i++)
			if (p[i]) {
				w=v[i];
				v[i]=g[i]*u[i]+f[i]*w;
				u[i]=f[i]*u[i]-g[i]*w;
				if (i < n) {
					v[i+1] -= v[i];
					u[i+1] -= u[i];
				}
			} else {
				aa=u[i+1];
				bb=v[i+1];
				u[i+1]=u[i]-(f[i]*aa-g[i]*bb);
				u[i]=aa;
				v[i+1]=v[i]-(g[i]*aa+f[i]*bb);
				v[i]=bb;
			}
		for (i=n; i>=1; i--) {
			i1=i+1;
			u[i]=(u[i]-matvec(i1,n,i,a,u)+
					(p[i] ? tamvec(i1,n,i,a,v) : a[i1][i]*v[i1]))/a[i][i];
			v[i]=(v[i]-matvec(i1,n,i,a,v)-
					(p[i] ? tamvec(i1,n,i,a,u) : a[i1][i]*u[i1]))/a[i][i];
		}
		w=1.0/sqrt(vecvec(1,n,0,u,u)+vecvec(1,n,0,v,v));
		for (j=1; j<=n; j++) {
			u[j] *= w;
			v[j] *= w;
		}
		count++;
	} while (w > tol);
	em[7]=w;
	em[9]=count;
	free_integer_vector(p,1);
	free_real_vector(g,1);
	free_real_vector(f,1);
}

