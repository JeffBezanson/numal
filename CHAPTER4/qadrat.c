#include "../real.h"


real_t qadrat(real_t *x, real_t a, real_t b, real_t (*fx)(real_t),
				real_t e[])
{
	real_t lint(real_t *, real_t (*)(real_t), real_t [], real_t, real_t,
				real_t, real_t, real_t, real_t, real_t, real_t, real_t,
				real_t, real_t, real_t, real_t, real_t);
	real_t f0,f2,f3,f5,f6,f7,f9,f14,hmin,hmax,re,ae,result;

	hmax=(b-a)/16.0;
	if (hmax == 0.0) return 0.0;
	re=e[1];
	ae=2.0*e[2]/fabs(b-a);
	e[3]=0.0;
	hmin=fabs(b-a)*re;
	*x=a;
	f0=(*fx)(*x);
	*x=a+hmax;
	f2=(*fx)(*x);
	*x=a+2.0*hmax;
	f3=(*fx)(*x);
	*x=a+4.0*hmax;
	f5=(*fx)(*x);
	*x=a+6.0*hmax;
	f6=(*fx)(*x);
	*x=a+8.0*hmax;
	f7=(*fx)(*x);
	*x=b-4.0*hmax;
	f9=(*fx)(*x);
	*x=b;
	f14=(*fx)(*x);
	result = lint(x,fx,e,a,b,f0,f2,f3,f5,f6,f7,f9,f14,
						hmin,hmax,re,ae)*16.0;
	return result;
}

real_t lint(real_t *x, real_t (*fx)(real_t), real_t e[], real_t x0,
				real_t xn, real_t f0, real_t f2, real_t f3, real_t f5,
				real_t f6, real_t f7, real_t f9, real_t f14,
				real_t hmin, real_t hmax, real_t re, real_t ae)
{
	/* this function is internally used by QADRAT */

	real_t v,w,h,xm,f1,f4,f8,f10,f11,f12,f13;

	xm=(x0+xn)/2.0;
	h=(xn-x0)/32.0;
	*x=xm+4.0*h;
	f8=(*fx)(*x);
	*x=xn-4.0*h;
	f11=(*fx)(*x);
	*x=xn-2.0*h;
	f12=(*fx)(*x);
	v=0.330580178199226*f7+0.173485115707338*(f6+f8)+
		0.321105426559972*(f5+f9)+0.135007708341042*(f3+f11)+
		0.165714514228223*(f2+f12)+0.393971460638127e-1*(f0+f14);
	*x=x0+h;
	f1=(*fx)(*x);
	*x=xn-h;
	f13=(*fx)(*x);
	w=0.260652434656970*f7+0.239063286684765*(f6+f8)+
		0.263062635477467*(f5+f9)+0.218681931383057*(f3+f11)+
		0.275789764664284e-1*(f2+f12)+0.105575010053846*(f1+f13)+
		0.157119426059518e-1*(f0+f14);
	if (fabs(h) < hmin) e[3] += 1.0;
	if (fabs(v-w) < fabs(w)*re+ae || fabs(h) < hmin)
		return h*w;
	else {
		*x=x0+6.0*h;
		f4=(*fx)(*x);
		*x=xn-6.0*h;
		f10=(*fx)(*x);
		v=0.245673430093324*f7+0.255786258286921*(f6+f8)+
			0.228526063690406*(f5+f9)+0.500557131525460e-1*(f4+f10)+
			0.177946487736780*(f3+f11)+0.584014599347449e-1*(f2+f12)+
			0.874830942871331e-1*(f1+f13)+
			0.189642078648079e-1*(f0+f14);
		return ((fabs(v-w) < fabs(v)*re+ae) ? h*v :
					(lint(x,fx,e,x0,xm,f0,f1,f2,f3,f4,f5,f6,f7,
							hmin,hmax,re,ae)-
					lint(x,fx,e,xn,xm,f14,f13,f12,f11,f10,f9,f8,f7,
							hmin,hmax,re,ae)));
	}
}

