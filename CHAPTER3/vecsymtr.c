#include "../real.h"


void vecsymtri(real_t d[], real_t b[], int n, int n1, int n2,
					real_t val[], real_t **vec, real_t em[])
{
	int *allocate_integer_vector(int, int);
	real_t *allocate_real_vector(int, int);
	void free_integer_vector(int *, int);
	void free_real_vector(real_t *, int);
	real_t vecvec(int, int, int, real_t [], real_t []);
	void elmveccol(int, int, int, real_t [], real_t **, real_t);
	real_t tamvec(int, int, int, real_t **, real_t []);
	int i,j,k,count,maxcount,countlim,orth,ind,iterate,*index;
	real_t bi,bi1,u,w,y,mi1,lambda,oldlambda,ortheps,valspread,spr,
			res,maxres,norm,newnorm,oldnorm,machtol,vectol,
			*m,*p,*q,*r,*x;

	index=allocate_integer_vector(1,n);
	m=allocate_real_vector(1,n);
	p=allocate_real_vector(1,n);
	q=allocate_real_vector(1,n);
	r=allocate_real_vector(1,n);
	x=allocate_real_vector(1,n);
	norm=em[1];
	machtol=em[0]*norm;
	valspread=em[4]*norm;
	vectol=em[6]*norm;
	countlim=em[8];
	ortheps=sqrt(em[0]);
	maxcount=ind=0;
	maxres=0.0;
	if (n1 > 1) {
		orth=em[5];
		oldlambda=val[n1-orth];
		for (k=n1-orth+1; k<=n1-1; k++) {
			lambda=val[k];
			spr=oldlambda-lambda;
			if (spr < machtol) lambda=oldlambda-machtol;
			oldlambda=lambda;
		}
	} else
		orth=1;
	for (k=n1; k<=n2; k++) {
		lambda=val[k];
		if (k > 1) {
			spr=oldlambda-lambda;
			if (spr < valspread) {
				if (spr < machtol) lambda=oldlambda-machtol;
				orth++;
			} else
				orth=1;
		}
		count=0;
		u=d[1]-lambda;
		bi=w=b[1];
		if (fabs(bi) < machtol) bi=machtol;
		for (i=1; i<=n-1; i++) {
			bi1=b[i+1];
			if (fabs(bi1) < machtol) bi1=machtol;
			if (fabs(bi) >= fabs(u)) {
				mi1=m[i+1]=u/bi;
				p[i]=bi;
				y=q[i]=d[i+1]-lambda;
				r[i]=bi1;
				u=w-mi1*y;
				w = -mi1*bi1;
				index[i]=1;
			} else {
				mi1=m[i+1]=bi/u;
				p[i]=u;
				q[i]=w;
				r[i]=0.0;
				u=d[i+1]-lambda-mi1*w;
				w=bi1;
				index[i]=0;
			}
			x[i]=1.0;
			bi=bi1;
		}	/* transform */
		p[n] = (fabs(u) < machtol) ? machtol : u;
		q[n]=r[n]=0.0;
		x[n]=1.0;
		iterate=1;
		while (iterate) {
			u=w=0.0;
			for (i=n; i>=1; i--) {
				y=u;
				u=x[i]=(x[i]-q[i]*u-r[i]*w)/p[i];
				w=y;
			}	/* next iteration */
			newnorm=sqrt(vecvec(1,n,0,x,x));
			if (orth > 1) {
				oldnorm=newnorm;
				for (j=k-orth+1; j<=k-1; j++)
					elmveccol(1,n,j,x,vec,-tamvec(1,n,j,vec,x));
				newnorm=sqrt(vecvec(1,n,0,x,x));
				if (newnorm < ortheps*oldnorm) {
					ind++;
					count=1;
					for (i=1; i<=ind-1; i++) x[i]=0.0;
					for (i=ind+1; i<=n; i++) x[i]=0.0;
					x[ind]=1.0;
					if (ind == n) ind=0;
					w=x[1];
					for (i=2; i<=n; i++) {
						if (index[i-1]) {
							u=w;
							w=x[i-1]=x[i];
						} else
							u=x[i];
						w=x[i]=u-m[i]*w;
					}
					continue;	/* iterate on */
				}	/* new start */
			}	/* orthogonalization */
			res=1.0/newnorm;
			if (res > vectol || count == 0) {
				count++;
				if (count <= countlim) {
					for (i=1; i<=n; i++) x[i] *= res;
					w=x[1];
					for (i=2; i<=n; i++) {
						if (index[i-1]) {
							u=w;
							w=x[i-1]=x[i];
						} else
							u=x[i];
						w=x[i]=u-m[i]*w;
					}
				} else
					break;
			} else
				break;
		}
		for (i=1; i<=n; i++) vec[i][k]=x[i]*res;
		if (count > maxcount) maxcount=count;
		if (res > maxres) maxres=res;
		oldlambda=lambda;
	}
	em[5]=orth;
	em[7]=maxres;
	em[9]=maxcount;
	free_integer_vector(index,1);
	free_real_vector(m,1);
	free_real_vector(p,1);
	free_real_vector(q,1);
	free_real_vector(r,1);
	free_real_vector(x,1);
}

