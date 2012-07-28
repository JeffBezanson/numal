#include "../real.h"
#include <stdlib.h>


real_t *allocate_real_vector(int, int);
real_t **allocate_real_matrix(int, int, int, int);
void free_real_vector(real_t *, int);
void free_real_matrix(real_t **, int, int, int);
void inivec(int, int, real_t [], real_t);
void inimat(int, int, int, int, real_t **, real_t);
void dupvec(int, int, int, real_t [], real_t []);
void dupmat(int, int, int, int, real_t **, real_t **);
void dupcolvec(int, int, int, real_t **, real_t []);
void mulrow(int, int, int, int, real_t **, real_t **, real_t);
void mulcol(int, int, int, int, real_t **, real_t **, real_t);
real_t vecvec(int, int, int, real_t [], real_t []);
real_t tammat(int, int, int, int, real_t **, real_t **);
real_t mattam(int, int, int, int, real_t **, real_t **);
void ichrowcol(int, int, int, int, real_t **);
void elmveccol(int, int, int, real_t [], real_t **, real_t);
int qrisngvaldec(real_t **, int, int, real_t [], real_t **, real_t []);
void praxismin(int, int, real_t *, real_t *, real_t *, int, int,
			   real_t [], real_t **, real_t *, real_t *, real_t *, real_t,	real_t,
			   real_t [], real_t [], int *, int *, real_t *, real_t, real_t,
			   real_t, real_t, real_t, real_t, real_t, real_t,
			   real_t (*)(int, real_t[]));

void praxis(int n, real_t x[], real_t (*funct)(int, real_t[]),
				real_t in[], real_t out[])
{
	int illc,i,j,k,k2,nl,maxf,nf,kl,kt,ktm,emergency;
	real_t s,sl,dn,dmin,fx,f1,lds,ldt,sf,df,qf1,qd0,qd1,qa,qb,qc,m2,m4,
			small,vsmall,large,vlarge,scbd,ldfac,t2,macheps,reltol,
			abstol,h,**v,*d,*y,*z,*q0,*q1,**a,em[8],l;

	d=allocate_real_vector(1,n);
	y=allocate_real_vector(1,n);
	z=allocate_real_vector(1,n);
	q0=allocate_real_vector(1,n);
	q1=allocate_real_vector(1,n);
	v=allocate_real_matrix(1,n,1,n);
	a=allocate_real_matrix(1,n,1,n);

	macheps=in[0];
	reltol=in[1];
	abstol=in[2];
	maxf=in[5];
	h=in[6];
	scbd=in[7];
	ktm=in[8];
	illc = in[9] < 0.0;
	small=macheps*macheps;
	vsmall=small*small;
	large=1.0/small;
	vlarge=1.0/vsmall;
	m2=reltol;
	m4=sqrt(m2);
	srand(1);
	ldfac = (illc ? 0.1 : 0.01);
	kt=nl=0;
	nf=1;
	out[3]=qf1=fx=(*funct)(n,x);
	abstol=t2=small+fabs(abstol);
	dmin=small;
	if (h < abstol*100.0) h=abstol*100;
	ldt=h;
	inimat(1,n,1,n,v,0.0);
	for (i=1; i<=n; i++) v[i][i]=1.0;
	d[1]=qd0=qd1=0.0;
	dupvec(1,n,0,q1,x);
	inivec(1,n,q0,0.0);
	emergency=0;

	while (1) {
		sf=d[1];
		d[1]=s=0.0;
		praxismin(1,2,&(d[1]),&s,&fx,0,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct);
		if (s <= 0.0) mulcol(1,n,1,1,v,v,-1.0);
		if (sf <= 0.9*d[1] || 0.9*sf >= d[1]) inivec(2,n,d,0.0);
		for (k=2; k<=n; k++) {
			dupvec(1,n,0,y,x);
			sf=fx;
			illc = (illc || kt > 0);
			while (1) {
				kl=k;
				df=0.0;
				if (illc) {
					/* random stop to get off resulting valley */
					for (i=1; i<=n; i++) {
						s=z[i]=(0.1*ldt+t2*pow(10.0,kt))*
									(rand()/(real_t)RAND_MAX-0.5);
						elmveccol(1,n,i,x,v,s);
					}
					fx=(*funct)(n,x);
					nf++;
				}
				for (k2=k; k2<=n; k2++) {
					sl=fx;
					s=0.0;
					praxismin(k2,2,&(d[k2]),&s,&fx,0,
						n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
						&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct);
					s = illc ? d[k2]*(s+z[k2])*(s+z[k2]) : sl-fx;
					if (df < s) {
						df=s;
						kl=k2;
					}
				}
				if (!illc && df < fabs(100.0*macheps*fx))
					illc=1;
				else
					break;
			}
			for (k2=1; k2<=k-1; k2++) {
				s=0.0;
				praxismin(k2,2,&(d[k2]),&s,&fx,0,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct);
			}
			f1=fx;
			fx=sf;
			lds=0.0;
			for (i=1; i<=n; i++) {
				sl=x[i];
				x[i]=y[i];
				y[i] = sl -= y[i];
				lds += sl*sl;
			}
			lds=sqrt(lds);
			if (lds > small) {
				for (i=kl-1; i>=k; i--) {
					for (j=1; j<=n; j++) v[j][i+1]=v[j][i];
					d[i+1]=d[i];
				}
				d[k]=0.0;
				dupcolvec(1,n,k,v,y);
				mulcol(1,n,k,k,v,v,1.0/lds);
				praxismin(k,4,&(d[k]),&lds,&f1,1,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct);
				if (lds <= 0.0) {
					lds = -lds;
					mulcol(1,n,k,k,v,v,-1.0);
				}
			}
			ldt *= ldfac;
			if (ldt < lds) ldt=lds;
			t2=m2*sqrt(vecvec(1,n,0,x,x))+abstol;
			kt = (ldt > 0.5*t2) ? 0 : kt+1;
			if (kt > ktm) {
				out[1]=0.0;
				emergency=1;
			}
		}
		if (emergency) break;
		/* quad */
		s=fx;
		fx=qf1;
		qf1=s;
		qd1=0.0;
		for (i=1; i<=n; i++) {
			s=x[i];
			x[i]=l=q1[i];
			q1[i]=s;
			qd1 += (s-l)*(s-l);
		}
		l=qd1=sqrt(qd1);
		s=0.0;
		if ((qd0*qd1 > FLT_MIN) && (nl >=3*n*n)) {
			praxismin(0,2,&s,&l,&qf1,1,
					n,x,v,&qa,&qb,&qc,qd0,qd1,q0,q1,&nf,
					&nl,&fx,m2,m4,dmin,ldt,reltol,abstol,small,h,funct);
			qa=l*(l-qd1)/(qd0*(qd0+qd1));
			qb=(l+qd0)*(qd1-l)/(qd0*qd1);
			qc=l*(l+qd0)/(qd1*(qd0+qd1));
		} else {
			fx=qf1;
			qa=qb=0.0;
			qc=1.0;
		}
		qd0=qd1;
		for (i=1; i<=n; i++) {
			s=q0[i];
			q0[i]=x[i];
			x[i]=qa*s+qb*x[i]+qc*q1[i];
		}
		/* end of quad */
		dn=0.0;
		for (i=1; i<=n; i++) {
			d[i]=1.0/sqrt(d[i]);
			if (dn < d[i]) dn=d[i];
		}
		for (j=1; j<=n; j++) {
			s=d[j]/dn;
			mulcol(1,n,j,j,v,v,s);
		}
		if (scbd > 1.0) {
			s=vlarge;
			for (i=1; i<=n; i++) {
				sl=z[i]=sqrt(mattam(1,n,i,i,v,v));
				if (sl < m4) z[i]=m4;
				if (s > sl) s=sl;
			}
			for (i=1; i<=n; i++) {
				sl=s/z[i];
				z[i]=1.0/sl;
				if (z[i] > scbd) {
					sl=1.0/scbd;
					z[i]=scbd;
				}
				mulrow(1,n,i,i,v,v,sl);
			}
		}
		for (i=1; i<=n; i++) ichrowcol(i+1,n,i,i,v);
		em[0]=em[2]=macheps;
		em[4]=10*n;
		em[6]=vsmall;
		dupmat(1,n,1,n,a,v);
		if (qrisngvaldec(a,n,n,d,v,em) != 0) {
			out[1]=2.0;
			emergency=1;
		}
		if (emergency) break;
		if (scbd > 1.0) {
			for (i=1; i<=n; i++) mulrow(1,n,i,i,v,v,z[i]);
			for (i=1; i<=n; i++) {
				s=sqrt(tammat(1,n,i,i,v,v));
				d[i] *= s;
				s=1.0/s;
				mulcol(1,n,i,i,v,v,s);
			}
		}
		for (i=1; i<=n; i++) {
			s=dn*d[i];
			d[i] = (s > large) ? vsmall :
						((s < small) ? vlarge : 1.0/(s*s));
		}
		/* sort */
		for (i=1; i<=n-1; i++) {
			k=i;
			s=d[i];
			for (j=i+1; j<=n; j++)
				if (d[j] > s) {
					k=j;
					s=d[j];
				}
			if (k > i) {
				d[k]=d[i];
				d[i]=s;
				for (j=1; j<=n; j++) {
					s=v[j][i];
					v[j][i]=v[j][k];
					v[j][k]=s;
				}
			}
		}
		/* end of sort */
		dmin=d[n];
		if (dmin < small) dmin=small;
		illc = (m2*d[1]) > dmin;
		if (nf >= maxf) {
			out[1]=1.0;
			break;
		}
	}
	out[2]=fx;
	out[4]=nf;
	out[5]=nl;
	out[6]=ldt;
	free_real_vector(d,1);
	free_real_vector(y,1);
	free_real_vector(z,1);
	free_real_vector(q0,1);
	free_real_vector(q1,1);
	free_real_matrix(v,1,n,1);
	free_real_matrix(a,1,n,1);
}

void praxismin(int j, int nits, real_t *d2, real_t *x1, real_t *f1,
			int fk, int n, real_t x[], real_t **v, real_t *qa, real_t *qb,
			real_t *qc, real_t qd0, real_t qd1, real_t q0[], real_t q1[],
			int *nf, int *nl, real_t *fx, real_t m2, real_t m4,
			real_t dmin, real_t ldt, real_t reltol, real_t abstol,
			real_t small, real_t h, real_t (*funct)(int, real_t[]))
{
	/* this function is internally used by PRAXIS */

	real_t praxisflin(real_t, int, int, real_t [], real_t **, real_t *,
					real_t *, real_t *, real_t, real_t, real_t [],
					real_t [], int *,	real_t (*)(int, real_t[]));
	int k,dz,loop;
	real_t x2,xm,f0,f2,fm,d1,t2,s,sf1,sx1;

	sf1 = *f1;
	sx1 = *x1;
	k=0;
	xm=0.0;
	f0 = fm = *fx;
	dz = *d2 < reltol;
	s=sqrt(vecvec(1,n,0,x,x));
	t2=m4*sqrt(fabs(*fx)/(dz ? dmin : *d2)+s*ldt)+m2*ldt;
	s=s*m4+abstol;
	if (dz && (t2 > s)) t2=s;
	if (t2 < small) t2=small;
	if (t2 > 0.01*h) t2=0.01*h;
	if (fk && (*f1 <= fm)) {
		xm = *x1;
		fm = *f1;
	}
	if (!fk || (fabs(*x1) < t2)) {
		*x1 = (*x1 > 0.0) ? t2 : -t2;
		*f1=praxisflin(*x1,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct);
	}
	if (*f1 <= fm) {
		xm = *x1;
		fm = *f1;
	}
	loop=1;
	while (loop) {
		if (dz) {
			/* evaluate praxisflin at another point and
				estimate the second derivative */
			x2 = (f0 < *f1) ? -(*x1) : (*x1)*2.0;
			f2=praxisflin(x2,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct);
			if (f2 <= fm) {
				xm=x2;
				fm=f2;
			}
			*d2=(x2*((*f1)-f0)-(*x1)*(f2-f0))/((*x1)*x2*((*x1)-x2));
		}
		/* estimate first derivative at 0 */
		d1=((*f1)-f0)/(*x1)-(*x1)*(*d2);
		dz=1;
		x2 = (*d2 <= small) ? ((d1 < 0.0) ? h : -h) : -0.5*d1/(*d2);
		if (fabs(x2) > h)	x2 = (x2 > 0.0) ? h : -h;
		while (1) {
			f2=praxisflin(x2,j,n,x,v,qa,qb,qc,qd0,qd1,q0,q1,nf,funct);
			if (k < nits && f2 > f0) {
				k++;
				if (f0 < *f1 && (*x1)*x2 > 0.0) break;
				x2=0.5*x2;
			} else {
				loop=0;
				break;
			}
		}
	}
	(*nl)++;
	if (f2 > fm)
		x2=xm;
	else
		fm=f2;
	*d2 = (fabs(x2*(x2-(*x1))) > small) ?
				((x2*((*f1)-f0)-(*x1)*(fm-f0))/((*x1)*x2*((*x1)-x2))) :
				((k > 0) ? 0.0 : *d2);
	if (*d2 <= small) *d2=small;
	*x1=x2;
	*fx=fm;
	if (sf1 < *fx) {
		*fx=sf1;
		*x1=sx1;
	}
	if (j > 0) elmveccol(1,n,j,x,v,*x1);
}

real_t praxisflin(real_t l, int j, int n, real_t x[], real_t **v,
			real_t *qa, real_t *qb, real_t *qc, real_t qd0, real_t qd1,
			real_t q0[], real_t q1[], int *nf,
			real_t (*funct)(int, real_t[]))
{
	/* this function is internally used by PRAXISMIN */

	int i;
	real_t *t,result;

	t=allocate_real_vector(1,n);
	if (j > 0)
		for (i=1; i<=n; i++) t[i]=x[i]+l*v[i][j];
	else {
		/* search along parabolic space curve */
		*qa=l*(l-qd1)/(qd0*(qd0+qd1));
		*qb=(l+qd0)*(qd1-l)/(qd0*qd1);
		*qc=l*(l+qd0)/(qd1*(qd0+qd1));
		for (i=1; i<=n; i++) t[i]=(*qa)*q0[i]+(*qb)*x[i]+(*qc)*q1[i];
	}
	(*nf)++;
	result=(*funct)(n,t);
	free_real_vector(t,1);
	return result;
}

