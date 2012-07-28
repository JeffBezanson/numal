#include "../real.h"


real_t minin(real_t *x, real_t *a, real_t *b, real_t (*fx)(real_t),
				real_t (*tolx)(real_t))
{
	real_t z,c,d,e,m,p,q,r,tol,t,u,v,w,fu,fv,fw,fz;

	c=(3.0-sqrt(5.0))/2.0;
	if (*a > *b) {
		z = *a;
		*a = *b;
		*b=z;
	}
	w = *x = *a;
	fw=(*fx)(*x);
	z = *x = *b;
	fz=(*fx)(*x);
	if (fz > fw) {
		z=w;
		w = *x;
		v=fz;
		fz=fw;
		fw=v;
	}
	v=w;
	fv=fw;
	e=0.0;
	while (1) {
		m=(*a + *b)*0.5;
		tol=(*tolx)(*x);
		t=tol*2.0;
		if (fabs(z-m) <= t-(*b - *a)*0.5) break;
		p=q=r=0.0;
		if (fabs(e) > tol) {
			r=(z-w)*(fz-fv);
			q=(z-v)*(fz-fw);
			p=(z-v)*q-(z-w)*r;
			q=(q-r)*2.0;
			if (q > 0.0)
				p = -p;
			else
				q = -q;
			r=e;
			e=d;
		}
		if (fabs(p) < fabs(q*r*0.5) && p > (*a-z)*q && p < (*b-z)*q) {
			d=p/q;
			u=z+d;
			if (u-(*a) < t || (*b)-u < t) d = ((z < m) ? tol : -tol);
		} else {
			e = ((z < m) ? *b : *a) - z;
			d=c*e;
		}
		u = *x = z + ((fabs(d) >= tol) ? d : ((d>0.0) ? tol : -tol));
		fu=(*fx)(*x);
		if (fu <= fz) {
			if (u < z)
				*b=z;
			else
				*a=z;
			v=w;
			fv=fw;
			w=z;
			fw=fz;
			z=u;
			fz=fu;
		} else {
			if (u < z)
				*a=u;
			else
				*b=u;
			if (fu <= fw) {
				v=w;
				fv=fw;
				w=u;
				fw=fu;
			} else
				if (fu <= fv || v == w) {
					v=u;
					fv=fu;
				}
		}
	}
	*x=z;
	return fz;
}
