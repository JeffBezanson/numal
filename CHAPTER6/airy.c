#include "../real.h"


void airy(real_t z, real_t *ai, real_t *aid, real_t *bi, real_t *bid,
			real_t *expon, int first)
{
	int n,l;
	real_t s,t,u,v,sc,tc,uc,vc,x,k1,k2,k3,k4,c,zt,si,co,expzt,sqrtz,wwl,
			pl,pl1,pl2,pl3;
	static real_t c1,c2,sqrt3,sqrt1opi,pio4,xx[11],ww[11];

	if (first) {
		sqrt3=1.73205080756887729;
		sqrt1opi=0.56418958354775629;
		pio4=0.78539816339744831;
		c1=0.355028053887817;
		c2=0.258819403792807;
		xx[1] =1.4083081072180964e1;
		xx[2] =1.0214885479197331e1;
		xx[3] =7.4416018450450930;
		xx[4] =5.3070943061781927;
		xx[5] =3.6340135029132462;
		xx[6] =2.3310652303052450;
		xx[7] =1.3447970842609268;
		xx[8] =6.4188858369567296e-1;
		xx[9] =2.0100345998121046e-1;
		xx[10]=8.0594359172052833e-3;
		ww[1] =3.1542515762964787e-14;
		ww[2] =6.6394210819584921e-11;
		ww[3] =1.7583889061345669e-8;
		ww[4] =1.3712392370435815e-6;
		ww[5] =4.4350966639284350e-5;
		ww[6] =7.1555010917718255e-4;
		ww[7] =6.4889566103335381e-3;
		ww[8] =3.6440415875773282e-2;
		ww[9] =1.4399792418590999e-1;
		ww[10]=8.1231141336261486e-1;
	}
	*expon=0.0;
	if (z >= -5.0 && z <= 8.0) {
		u=v=t=uc=vc=tc=1.0;
		s=sc=0.5;
		n=3;
		x=z*z*z;
		while (fabs(u)+fabs(v)+fabs(s)+fabs(t) > 1.0e-18) {
			u=u*x/(n*(n-1));
			v=v*x/(n*(n+1));
			s=s*x/(n*(n+2));
			t=t*x/(n*(n-2));
			uc += u;
			vc += v;
			sc += s;
			tc += t;
			n += 3;
		}
		*bi=sqrt3*(c1*uc+c2*z*vc);
		*bid=sqrt3*(c1*z*z*sc+c2*tc);
		if (z < 2.5) {
			*ai=c1*uc-c2*z*vc;
			*aid=c1*sc*z*z-c2*tc;
			return;
		}
	}
	k1=k2=k3=k4=0.0;
	sqrtz=sqrt(fabs(z));
	zt=0.666666666666667*fabs(z)*sqrtz;
	c=sqrt1opi/sqrt(sqrtz);
	if (z < 0.0) {
		z = -z;
		co=cos(zt-pio4);
		si=sin(zt-pio4);
		for (l=1; l<=10; l++) {
			wwl=ww[l];
			pl=xx[l]/zt;
			pl2=pl*pl;
			pl1=1.0+pl2;
			pl3=pl1*pl1;
			k1 += wwl/pl1;
			k2 += wwl*pl/pl1;
			k3 += wwl*pl*(1.0+pl*(2.0/zt+pl))/pl3;
			k4 += wwl*(-1.0-pl*(1.0+pl*(zt-pl))/zt)/pl3;
		}
		*ai=c*(co*k1+si*k2);
		*aid=0.25*(*ai)/z-c*sqrtz*(co*k3+si*k4);
		*bi=c*(co*k2-si*k1);
		*bid=0.25*(*bi)/z-c*sqrtz*(co*k4-si*k3);
	} else {
		if (z < 9.0)
			expzt=exp(zt);
		else {
			expzt=1.0;
			*expon=zt;
		}
		for (l=1; l<=10; l++) {
			wwl=ww[l];
			pl=xx[l]/zt;
			pl1=1.0+pl;
			pl2=1.0-pl;
			k1 += wwl/pl1;
			k2 += wwl*pl/(zt*pl1*pl1);
			k3 += wwl/pl2;
			k4 += wwl*pl/(zt*pl2*pl2);
		}
		*ai=0.5*c*k1/expzt;
		*aid=(*ai)*(-0.25/z-sqrtz)+0.5*c*sqrtz*k2/expzt;
		if (z >= 8.0) {
			*bi=c*k3*expzt;
			*bid=(*bi)*(sqrtz-0.25/z)-c*k4*sqrtz*expzt;
		}
	}
}
