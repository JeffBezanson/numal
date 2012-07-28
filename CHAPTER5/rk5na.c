#include "../real.h"


void rk5na(real_t x[], real_t xa[], real_t (*b)(int, real_t[]),
			real_t (*fxj)(int, int, real_t[]), real_t e[], real_t d[],
			int fi, int n, int l, int pos)
{
	real_t *allocate_real_vector(int, int);
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_vector(real_t *, int);
	void free_real_matrix(real_t **, int, int, int);
	void rk5narkstep(real_t, int, int, real_t,
					real_t (*)(int, int, real_t[]), real_t [], real_t [],
					real_t [], real_t [], real_t **);
	int j,i,first,fir,rej,t,ext,extrapolate;
	real_t fhm,s,s0,cond0,s1,cond1,h,absh,tol,fh,hl,mu,mu1,
			fzero,*y,*xl,*discr,**k,e1[3],
			c,fc,bb,fb,a,fa,dd,fd,fdb,fda,w,mb,m,p,q;

	y=allocate_real_vector(0,n);
	xl=allocate_real_vector(0,n);
	discr=allocate_real_vector(0,n);
	k=allocate_real_matrix(0,5,0,n);
	if (fi) {
		for (i=0; i<=n; i++) d[i+3]=xa[i];
		d[1]=d[2]=0.0;
	}
	for (i=0; i<=n; i++) x[i]=xl[i]=d[i+3];
	s=d[1];
	first=fir=1;
	h=e[0]+e[1];
	for (i=1; i<=n; i++) {
		absh=e[2*i]+e[2*i+1];
		if (h > absh) h=absh;
	}
	if (fi) {
		j=l;
		t=(*fxj)(n,j,x)*h < 0.0;
		if ((t && pos) || !(t || pos)) h = -h;
	} else
		if (d[2]*h < 0.0) h = -h;
	i=0;
	while (1) {
		rk5narkstep(h,i,n,mu,fxj,x,xl,y,discr,k);
		rej=0;
		fhm=0.0;
		absh=fabs(h);
		for (i=0; i<=n; i++) {
			tol=e[2*i]*fabs(k[0][i])+e[2*i+1]*absh;
			rej=(tol < discr[i] || rej);
			fh=discr[i]/tol;
			if (fh > fhm) fhm=fh;
		}
		mu=1.0/(1.0+fhm)+0.45;
		if (rej) {
			h *= mu;
			i=1;
		} else {
			if (first) {
				first=fir;
				hl=h;
				h *= mu;
			} else {
				fh=mu*h/hl+mu-mu1;
				hl=h;
				h *= fh;
			}
			rk5narkstep(hl,2,n,mu,fxj,x,xl,y,discr,k);
			mu1=mu;
			s += hl;
			if (fir) {
				cond0=(*b)(n,x);
				fir=0;
				if (!fi) h=d[2];
			} else {
				d[2]=h;
				cond1=(*b)(n,x);
				if (cond0*cond1 <= 0.0) break;
				cond0=cond1;
			}
			for (i=0; i<=n; i++) d[i+3]=xl[i]=x[i];
			d[1]=s0=s;
			i=0;
		}
	}
	e1[1]=e[2*n+2];
	e1[2]=e[2*n+3];
	s1=s;
	s=s0;
	/* find zero */
	bb=s;
	if (s == s0)
		fzero=cond0;
	else
		if (s == s1)
			fzero=cond1;
		else {
			rk5narkstep(s-s0,3,n,mu,fxj,x,xl,y,discr,k);
			fzero=(*b)(n,x);
		}
	fb=fzero;
	a=s=s1;
	if (s == s0)
		fzero=cond0;
	else
		if (s == s1)
			fzero=cond1;
		else {
			rk5narkstep(s-s0,3,n,mu,fxj,x,xl,y,discr,k);
			fzero=(*b)(n,x);
		}
	fa=fzero;
	c=a;
	fc=fa;
	ext=0;
	extrapolate=1;
	while (extrapolate) {
		if (fabs(fc) < fabs(fb)) {
			if (c != a) {
				dd=a;
				fd=fa;
			}
			a=bb;
			fa=fb;
			bb=s=c;
			fb=fc;
			c=a;
			fc=fa;
		}
		tol=fabs(e1[1]*s)+fabs(e1[2]);
		m=(c+bb)*0.5;
		mb=m-bb;
		if (fabs(mb) > tol) {
			if (ext > 2)
				w=mb;
			else {
				if (mb == 0.0)
					tol=0.0;
				else
					if (mb < 0.0) tol = -tol;
				p=(bb-a)*fb;
				if (ext <= 1)
					q=fa-fb;
				else {
					fdb=(fd-fb)/(dd-bb);
					fda=(fd-fa)/(dd-a);
					p *= fda;
					q=fdb*fa-fda*fb;
				}
				if (p < 0.0) {
					p = -p;
					q = -q;
				}
				w=(p<FLT_MIN || p<=q*tol) ? tol : ((p<mb*q) ? p/q : mb);
			}
			dd=a;
			fd=fa;
			a=bb;
			fa=fb;
			s = bb += w;
			if (s == s0)
				fzero=cond0;
			else
				if (s == s1)
					fzero=cond1;
				else {
					rk5narkstep(s-s0,3,n,mu,fxj,x,xl,y,discr,k);
					fzero=(*b)(n,x);
				}
			fb=fzero;
			if ((fc >= 0.0) ? (fb >= 0.0) : (fb <= 0.0)) {
				c=a;
				fc=fa;
				ext=0;
			} else
				ext = (w == mb) ? 0 : ext+1;
		} else
			break;
	}
	/* end of finding zero */
	rk5narkstep(s-s0,3,n,mu,fxj,x,xl,y,discr,k);
	for (i=0; i<=n; i++) d[i+3]=x[i];
	d[1]=s;
	free_real_vector(y,0);
	free_real_vector(xl,0);
	free_real_vector(discr,0);
	free_real_matrix(k,0,5,0);
}

void rk5narkstep(real_t h, int d, int n, real_t mu,
				real_t (*fxj)(int, int, real_t[]), real_t x[], real_t xl[],
				real_t y[], real_t discr[], real_t **k)
{
	/* this function is internally used by RK5NA */

	int i,j;
	real_t p,s;

	if (d != 2) {
		if (d == 1)
			for (i=0; i<=n; i++) k[0][i] *= mu;
		else {
			for (i=0; i<=n; i++) x[i]=xl[i];
			for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
			s=0.0;
			for (j=0; j<=n; j++) s += y[j]*y[j];
			p=h/sqrt(s);
			for (i=0; i<=n; i++) k[0][i]=y[i]*p;
		}
		for (i=0; i<=n; i++)	x[i]=xl[i]+k[0][i]/4.5;
		for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
		s=0.0;
		for (j=0; j<=n; j++) s += y[j]*y[j];
		p=h/sqrt(s);
		for (i=0; i<=n; i++) k[1][i]=y[i]*p;
		for (i=0; i<=n; i++)	x[i]=xl[i]+(k[0][i]+k[1][i]*3.0)/12.0;
		for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
		s=0.0;
		for (j=0; j<=n; j++) s += y[j]*y[j];
		p=h/sqrt(s);
		for (i=0; i<=n; i++) k[2][i]=y[i]*p;
		for (i=0; i<=n; i++)	x[i]=xl[i]+(k[0][i]+k[2][i]*3.0)/8.0;
		for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
		s=0.0;
		for (j=0; j<=n; j++) s += y[j]*y[j];
		p=h/sqrt(s);
		for (i=0; i<=n; i++) k[3][i]=y[i]*p;
		for (i=0; i<=n; i++)
			x[i]=xl[i]+(k[0][i]*53.0-k[1][i]*135.0+k[2][i]*126.0+
							k[3][i]*56.0)/125.0;
		for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
		s=0.0;
		for (j=0; j<=n; j++) s += y[j]*y[j];
		p=h/sqrt(s);
		for (i=0; i<=n; i++) k[4][i]=y[i]*p;
		if (d <= 1) {
			for (i=0; i<=n; i++)
				x[i]=xl[i]+(k[0][i]*133.0-k[1][i]*378.0+k[2][i]*276.0+
								k[3][i]*112.0+k[4][i]*25.0)/168.0;
			for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
			s=0.0;
			for (j=0; j<=n; j++) s += y[j]*y[j];
			p=h/sqrt(s);
			for (i=0; i<=n; i++) k[5][i]=y[i]*p;
			for (i=0; i<=n; i++)
				discr[i]=fabs(k[0][i]*21.0-k[2][i]*162.0+
						k[3][i]*224.0-k[4][i]*125.0+k[5][i]*42.0)/14.0;
			return;
		}
	}
	for (i=0; i<=n; i++)
		x[i]=xl[i]+(-k[0][i]*63.0+k[1][i]*189.0-k[2][i]*36.0-
						k[3][i]*112.0+k[4][i]*50.0)/28.0;
	for (j=0; j<=n; j++) y[j]=(*fxj)(n,j,x);
	s=0.0;
	for (j=0; j<=n; j++) s += y[j]*y[j];
	p=h/sqrt(s);
	for (i=0; i<=n; i++) k[5][i]=y[i]*p;
	for (i=0; i<=n; i++)
		x[i]=xl[i]+(k[0][i]*35.0+k[2][i]*162.0+k[4][i]*125.0+
						k[5][i]*14.0)/336.0;
}
