#include "../real.h"


void rk4a(real_t *x, real_t xa, real_t (*b)(real_t, real_t),
			real_t *y, real_t ya, real_t (*fxy)(real_t, real_t),
			real_t e[], real_t d[], int fi, int xdir, int pos)
{
	void rk4arkstep(real_t *, real_t, real_t, real_t *, real_t, real_t,
				real_t (*)(real_t, real_t), int, int, real_t *, real_t *,
				real_t *, real_t *, real_t *,	real_t *, real_t *, real_t);
	int i,iv,first,fir,rej,t,next,ext,extrapolate;
	real_t k0,k1,k2,k3,k4,k5,fhm,absh,discr,s,xl,cond0,s1,cond1,
			yl,hmin,h,zl,tol,hl,mu,mu1,e1[3],fzero,
			c,fc,bb,fb,a,fa,dd,fd,fdb,fda,w,mb,m,p,q;

	if (fi) {
		d[3]=xa;
		d[4]=ya;
		d[0]=1.0;
	}
	d[1]=0.0;
	*x=xl=d[3];
	*y=yl=d[4];
	iv = d[0] > 0.0;
	first=fir=1;
	hmin=e[0]+e[1];
	h=e[2]+e[3];
	if (h < hmin) hmin=h;
	while (1) {
		zl=(*fxy)(*x,*y);
		if (fabs(zl) <= 1.0) {
			if (!iv) {
				d[2] = h /= zl;
				d[0]=1.0;
				iv=first=1;
			}
			if (fir) {
				t=(((iv && xdir) || (!(iv || xdir))) ? h : h*zl) < 0.0;
				if (fi ? ((t && pos) || (!(t || pos))) : (h*d[2] < 0))
							h = -h;
			}
			i=1;
		} else {
			if (iv) {
				if (!fir) d[2] = h *= zl;
				d[0] = -1.0;
				iv=0;
				first=1;
			}
			if (fir) {
				h=e[0]+e[1];
				t=(((iv && xdir) || (!(iv || xdir))) ? h : h*zl) < 0.0;
				if (fi ? ((t && pos) || (!(t || pos))) : (h*d[2] < 0))
							h = -h;
			}
			i=1;
		}
		next=0;
		while (1) {
			absh=fabs(h);
			if (absh < hmin) {
				h = (h == 0.0) ? 0.0 : ((h > 0.0) ? hmin : -hmin);
				absh=hmin;
			}
			if (iv) {
				rk4arkstep(x,xl,h,y,yl,zl,fxy,i,0,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				tol=e[2]*fabs(k0)+e[3]*absh;
			} else {
				rk4arkstep(y,yl,h,x,xl,1.0/zl,fxy,i,1,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				tol=e[0]*fabs(k0)+e[1]*absh;
			}
			rej = discr > tol;
			mu=tol/(tol+discr)+0.45;
			if (!rej) break;
			if (absh <= hmin) {
				if (iv) {
					*x=xl+h;
					*y=yl+k0;
				} else {
					*x=xl+k0;
					*y=yl+h;
				}
				d[1] += 1.0;
				first=1;
				next=1;
				break;
			}
			h *= mu;
			i=0;
		}
		if (!next) {
			if (first) {
				first=fir;
				hl=h;
				h *= mu;
			} else {
				fhm=mu*h/hl+mu-mu1;
				hl=h;
				h *= fhm;
			}
			if (iv)
				rk4arkstep(x,xl,hl,y,yl,zl,fxy,2,0,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
			else
				rk4arkstep(y,yl,hl,x,xl,zl,fxy,2,1,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
			mu1=mu;
		}
		if (fir) {
			fir=0;
			cond0=(*b)(*x,*y);
			if (!(fi || rej)) h=d[2];
		} else {
			d[2]=h;
			cond1=(*b)(*x,*y);
			if (cond0*cond1 <= 0.0) break;
			cond0=cond1;
		}
		d[3]=xl=(*x);
		d[4]=yl=(*y);
	}
	e1[1]=e[4];
	e1[2]=e[5];
	s1 = iv ? *x : *y;
	s = iv ? xl : yl;
	/* find zero */
	bb=s;
	if (iv) {
		if (s == xl)
			fzero=cond0;
		else
			if (s == s1)
				fzero=cond1;
			else {
				rk4arkstep(x,xl,s-xl,y,yl,zl,fxy,3,0,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				fzero=(*b)(*x,*y);
			}
	} else {
		if (s == yl)
			fzero=cond0;
		else
			if (s == s1)
				fzero=cond1;
			else {
				rk4arkstep(y,yl,s-yl,x,xl,zl,fxy,3,1,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				fzero=(*b)(*x,*y);
			}
	}
	fb=fzero;
	a=s=s1;
	if (iv) {
		if (s == xl)
			fzero=cond0;
		else
			if (s == s1)
				fzero=cond1;
			else {
				rk4arkstep(x,xl,s-xl,y,yl,zl,fxy,3,0,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				fzero=(*b)(*x,*y);
			}
	} else {
		if (s == yl)
			fzero=cond0;
		else
			if (s == s1)
				fzero=cond1;
			else {
				rk4arkstep(y,yl,s-yl,x,xl,zl,fxy,3,1,
						&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
				fzero=(*b)(*x,*y);
			}
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
			if (iv) {
				if (s == xl)
					fzero=cond0;
				else
					if (s == s1)
						fzero=cond1;
					else {
						rk4arkstep(x,xl,s-xl,y,yl,zl,fxy,3,0,
								&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
						fzero=(*b)(*x,*y);
					}
			} else {
				if (s == yl)
					fzero=cond0;
				else
					if (s == s1)
						fzero=cond1;
					else {
						rk4arkstep(y,yl,s-yl,x,xl,zl,fxy,3,1,
								&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
						fzero=(*b)(*x,*y);
					}
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
	s1 = iv ? *x : *y;
	if (iv)
		rk4arkstep(x,xl,s-xl,y,yl,zl,fxy,3,0,
				&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
	else
		rk4arkstep(y,yl,s-yl,x,xl,zl,fxy,3,1,
				&k0,&k1,&k2,&k3,&k4,&k5,&discr,mu);
	d[3]=(*x);
	d[4]=(*y);
}

void rk4arkstep(real_t *x, real_t xl, real_t h, real_t *y, real_t yl,
			real_t zl, real_t (*fxy)(real_t, real_t), int d,	int invf,
			real_t *k0, real_t *k1, real_t *k2, real_t *k3, real_t *k4,
			real_t *k5, real_t *discr, real_t mu)
{
	/* this function is internally used by RK4A */

	if (d != 2) {
		if (d == 3) {
			*x=xl;
			*y=yl;
			*k0=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
		} else
			if (d == 1)
				*k0=zl*h;
			else
				*k0 *= mu;
		*x=xl+h/4.5;
		*y=yl+(*k0)/4.5;
		*k1=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
		*x=xl+h/3.0;
		*y=yl+((*k0)+(*k1)*3.0)/12.0;
		*k2=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
		*x=xl+h*0.5;
		*y=yl+((*k0)+(*k2)*3.0)/8.0;
		*k3=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
		*x=xl+h*0.8;
		*y=yl+((*k0)*53.0-(*k1)*135.0+(*k2)*126.0+(*k3)*56.0)/125.0;
		*k4=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
		if (d <= 1) {
			*x=xl+h;
			*y=yl+((*k0)*133.0-(*k1)*378.0+(*k2)*276.0+
					(*k3)*112.0+(*k4)*25.0)/168.0;
			*k5=(invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
			*discr=fabs((*k0)*21.0-(*k2)*162.0+(*k3)*224.0-
						(*k4)*125.0+(*k5)*42.0)/14.0;
			return;
		}
	}
	*x=xl+h;
	*y=yl+(-(*k0)*63.0+(*k1)*189.0-(*k2)*36.0-
			(*k3)*112.0+(*k4)*50.0)/28.0;
	*k5 = (invf ? (1.0/(*fxy)(*x,*y)) : (*fxy)(*x,*y))*h;
	*y=yl+((*k0)*35.0+(*k2)*162.0+(*k4)*125.0+(*k5)*14.0)/336.0;
}
