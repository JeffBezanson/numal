#include "../real.h"


void nonexpbessk01(real_t x, real_t *k0, real_t *k1)
{
	if (x <= 1.5) {
		void bessk01(real_t, real_t *, real_t *);
		real_t expx;
		expx=exp(x);
		bessk01(x,k0,k1);
		*k0 *= expx;
		*k1 *= expx;
	} else if (x <= 5.0) {
		int i,r;
		real_t t2,s1,s2,term1,term2,sqrtexpr,exph2,x2;
		static real_t fac[20]={0.90483741803596, 0.67032004603564,
			0.40656965974060, 0.20189651799466, 0.82084998623899e-1,
			0.27323722447293e-1, 0.74465830709243e-2,
			0.16615572731739e-2, 0.30353913807887e-3,
			0.45399929762485e-4, 0.55595132416500e-5,
			0.55739036926944e-6, 0.45753387694459e-7,
			0.30748798795865e-8, 0.16918979226151e-9,
			0.76218651945127e-11, 0.28111852987891e-12,
			0.84890440338729e-14, 0.2098791048793e-15,
			0.42483542552916e-17};
		s1=0.5;
		s2=0.0;
		r=0.0;
		x2=x+x;
		exph2=1.0/sqrt(5.0*x);
		for (i=0; i<=19; i++) {
			r += 1.0;
			t2=r*r/10.0;
			sqrtexpr=sqrt(t2/x2+1.0);
			term1=fac[i]/sqrtexpr;
			term2=fac[i]*sqrtexpr*t2;
			s1 += term1;
			s2 += term2;
		}
		*k0 = exph2*s1;
		*k1 = exph2*s2*2.0;
	} else {
		int r,i;
		real_t br,br1,br2,cr,cr1,cr2,ermin1,erplus1,er,f0,f1,
				expx,y,y2;
		static real_t dr[14]={0.27545e-15, -0.172697e-14,
				0.1136042e-13, -0.7883236e-13, 0.58081063e-12,
				-0.457993633e-11, 0.3904375576e-10, -0.36454717921e-9,
				0.379299645568e-8, -0.450473376411e-7,
				0.63257510850049e-6, -0.11106685196665e-4,
				0.26953261276272e-3, -0.11310504646928e-1};
		y=10.0/x-1.0;
		y2=y+y;
		r=30;
		br1=br2=cr1=cr2=erplus1=er=0.0;
		for (i=0; i<=13; i++) {
			r -= 2;
			br=y2*br1-br2+dr[i];
			cr=cr1*y2-cr2+er;
			ermin1=r*dr[i]+erplus1;
			erplus1=er;
			er=ermin1;
			br2=br1;
			br1=br;
			cr2=cr1;
			cr1=cr;
		}
		f0=y*br1-br2+0.9884081742308258;
		f1=y*cr1-cr2+er/2.0;
		expx=sqrt(1.5707963267949/x);
		*k0 = f0 *=expx;
		*k1 = (1.0+0.5/x)*f0+(10.0/x/x)*expx*f1;
	}
}
