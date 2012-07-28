#include "../real.h"


real_t integral(real_t a, real_t b, real_t (*fx)(real_t), real_t e[],
					int ua, int ub)
{
	real_t integralqad(int, real_t (*)(real_t), real_t [], real_t *,
				real_t *, real_t *, real_t *, real_t *, real_t *, real_t,
				real_t, real_t);
	real_t x0,x1,x2,f0,f1,f2,re,ae,b1,x;

	re=e[1];
	if (ub)
		ae=e[2]*180.0/fabs(b-a);
	else
		ae=e[2]*90.0/fabs(b-a);
	if (ua) {
		e[3]=e[4]=0.0;
		x=x0=a;
		f0=(*fx)(x);
	} else {
		x=x0=a=e[5];
		f0=e[6];
	}
	e[5]=x=x2=b;
	e[6]=f2=(*fx)(x);
	e[4] += integralqad(0,fx,e,&x0,&x1,&x2,&f0,&f1,&f2,re,ae,b1);
	if (!ub) {
		if (a < b) {
			b1=b-1.0;
			x0=1.0;
		} else {
			b1=b+1.0;
			x0 = -1.0;
		}
		f0=e[6];
		e[5]=x2=0.0;
		e[6]=f2=0.0;
		ae=e[2]*90.0;
		e[4] -= integralqad(1,fx,e,&x0,&x1,&x2,&f0,&f1,&f2,re,ae,b1);
	}
	return e[4];
}

real_t integralqad(int transf, real_t (*fx)(real_t), real_t e[],
			real_t *x0, real_t *x1, real_t *x2, real_t *f0, real_t *f1,
			real_t *f2, real_t re, real_t ae, real_t b1)
{
	/* this function is internally used by INTEGRAL */

	void integralint(int, real_t (*)(real_t), real_t [], real_t *,
			real_t *, real_t *, real_t *, real_t *, real_t *,
			real_t *, real_t, real_t, real_t, real_t);
	real_t sum,hmin,x,z;

	hmin=fabs((*x0)-(*x2))*re;
	x=(*x1)=((*x0)+(*x2))*0.5;
	if (transf) {
		z=1.0/x;
		x=z+b1;
		(*f1)=(*fx)(x)*z*z;
	} else
		(*f1)=(*fx)(x);
	sum=0.0;
	integralint(transf,fx,e,x0,x1,x2,f0,f1,f2,&sum,re,ae,b1,hmin);
	return sum/180.0;
}

void integralint(int transf, real_t (*fx)(real_t), real_t e[],
			real_t *x0, real_t *x1, real_t *x2, real_t *f0, real_t *f1,
			real_t *f2, real_t *sum, real_t re, real_t ae, real_t b1,
			real_t hmin)
{
	/* this function is internally used by INTEGRALQAD of INTEGRAL */

	int anew;
	real_t x3,x4,f3,f4,h,x,z,v,t;

	x4=(*x2);
	(*x2)=(*x1);
	f4=(*f2);
	(*f2)=(*f1);
	anew=1;
	while (anew) {
		anew=0;
		x=(*x1)=((*x0)+(*x2))*0.5;
		if (transf) {
			z=1.0/x;
			x=z+b1;
			(*f1)=(*fx)(x)*z*z;
		} else
			(*f1)=(*fx)(x);
		x=x3=((*x2)+x4)*0.5;
		if (transf) {
			z=1.0/x;
			x=z+b1;
			f3=(*fx)(x)*z*z;
		} else
			f3=(*fx)(x);
		h=x4-(*x0);
		v=(4.0*((*f1)+f3)+2.0*(*f2)+(*f0)+f4)*15.0;
		t=6.0*(*f2)-4.0*((*f1)+f3)+(*f0)+f4;
		if (fabs(t) < fabs(v)*re+ae)
			(*sum) += (v-t)*h;
		else if (fabs(h) < hmin)
			e[3] += 1.0;
		else {
			integralint(transf,fx,e,x0,x1,x2,f0,f1,f2,sum,
							re,ae,b1,hmin);
			*x2=x3;
			*f2=f3;
			anew=1;
		}
		if (!anew) {
			*x0=x4;
			*f0=f4;
		}
	}
}

