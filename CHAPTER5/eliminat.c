#include "../real.h"


void elimination(real_t **u, int lj, int uj, int ll, int ul,
			void (*residual)(int, int, int, int, real_t **),
			real_t a, real_t b, int *n, real_t discr[], int *k,
			real_t *rateconv, real_t *domeigval,
			void (*out)(real_t **, int, int, int, int, int *, real_t [],
							int, real_t, real_t))
{
	real_t **allocate_real_matrix(int, int, int, int);
	void free_real_matrix(real_t **, int, int, int);
	void richardson(real_t **, int, int, int, int,
			int, void (*)(int, int, int, int, real_t **),
			real_t, real_t, int *, real_t [], int *, real_t *, real_t *,
			void (*)(real_t **, int, int, int, int, int *, real_t [],
							int, real_t, real_t));
	real_t optpol(real_t, real_t, real_t, real_t);
	int ext,extrapolate;
	real_t pi,auxcos,c,d,
			cc,fc,bb,fb,aa,fa,dd,fd,fdb,fda,w,mb,tol,m,p,q;

	pi=3.14159265358979;
	c=1.0;
	if (optpol(c,a,b,*domeigval) < 0.0) {
		d=0.5*pi*sqrt(fabs(b/(*domeigval)));
		while (1) {
			d += d;
			/* finding zero */
			bb=c;
			fb=optpol(c,a,b,*domeigval);
			aa=c=d;
			fa=optpol(c,a,b,*domeigval);
			cc=aa;
			fc=fa;
			ext=0;
			extrapolate=1;
			while (extrapolate) {
				if (fabs(fc) < fabs(fb)) {
					if (cc != aa) {
						dd=aa;
						fd=fa;
					}
					aa=bb;
					fa=fb;
					bb=c=cc;
					fb=fc;
					cc=aa;
					fc=fa;
				}
				tol=c*1.0e-3;
				m=(cc+bb)*0.5;
				mb=m-bb;
				if (fabs(mb) > tol) {
					if (ext > 2)
						w=mb;
					else {
						if (mb == 0.0)
							tol=0.0;
						else
							if (mb < 0.0) tol = -tol;
						p=(bb-aa)*fb;
						if (ext <= 1)
							q=fa-fb;
						else {
							fdb=(fd-fb)/(dd-bb);
							fda=(fd-fa)/(dd-aa);
							p *= fda;
							q=fdb*fa-fda*fb;
						}
						if (p < 0.0) {
							p = -p;
							q = -q;
						}
						w=(p<FLT_MIN || p<=q*tol) ? tol :
								((p<mb*q) ? p/q : mb);
					}
					dd=aa;
					fd=fa;
					aa=bb;
					fa=fb;
					c = bb += w;
					fb=optpol(c,a,b,*domeigval);
					if ((fc >= 0.0) ? (fb >= 0.0) : (fb <= 0.0)) {
						cc=aa;
						fc=fa;
						ext=0;
					} else
						ext = (w == mb) ? 0 : ext+1;
				} else
					break;
			}
			d=cc;
			if ((fc >= 0.0) ? (fb <= 0.0) : (fb >= 0.0)) {
				*n=floor(c+0.5);
				break;
			}
		}
	} else
		*n=1;
	auxcos=cos(0.5*pi/(*n));
	richardson(u,lj,uj,ll,ul,1,residual,
		(2.0*(*domeigval)+b*(auxcos-1.0))/(auxcos+1.0),b,n,discr,k,
		rateconv,domeigval,out);
}

real_t optpol(real_t x, real_t a, real_t b, real_t domeigval)
{
	/* this function is internally used by ELIMINATION */

	real_t pi,w,y;

	pi=3.14159265358979;
	w=(b*cos(0.5*pi/x)+domeigval)/(b-domeigval);
	if (w < -1.0) w = -1.0;
	if (fabs(w) <= 1.0) {
		y=acos(w);
		return 2.0*sqrt(a/b)+tan(x*y)*(y-b*pi*sin(0.5*pi/x)*
					0.5/(x*(b-domeigval)*sqrt(fabs(1.0-w*w))));
	} else {
		y=log(w+sqrt(fabs(w*w-1.0)));
		return 2.0*sqrt(a/b)-tanh(x*y)*(y+b*pi*sin(0.5*pi/x)*
					0.5/(x*(b-domeigval)*sqrt(fabs(w*w-1.0))));
	}
}
