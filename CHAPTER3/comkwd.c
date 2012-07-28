#include "../real.h"


void comkwd(real_t pr, real_t pi, real_t qr, real_t qi,
				real_t *gr, real_t *gi, real_t *kr, real_t *ki)
{
	void commul(real_t, real_t, real_t, real_t, real_t *, real_t *);
	void comdiv(real_t, real_t, real_t, real_t, real_t *, real_t *);
	void comsqrt(real_t, real_t, real_t *, real_t *);
	real_t hr,hi;

	if (qr == 0.0 && qi == 0.0) {
		*kr = *ki = 0.0;
		*gr = pr*2.0;
		*gi = pi*2.0;
		return;
	}
	if (pr == 0.0 && pi == 0.0) {
		comsqrt(qr,qi,gr,gi);
		*kr = -(*gr);
		*ki = -(*gi);
		return;
	}
	if (fabs(pr) > 1.0 || fabs(pi) > 1.0) {
		comdiv(qr,qi,pr,pi,&hr,&hi);
		comdiv(hr,hi,pr,pi,&hr,&hi);
		comsqrt(1.0+hr,hi,&hr,&hi);
		commul(pr,pi,hr+1.0,hi,gr,gi);
	} else {
		comsqrt(qr+(pr+pi)*(pr-pi),qi+pr*pi*2.0,&hr,&hi);
		if (pr*hr+pi*hi > 0.0) {
			*gr = pr+hr;
			*gi = pi+hi;
		} else {
			*gr = pr-hr;
			*gi = pi-hi;
		}
	}
	comdiv(-qr,-qi,*gr,*gi,kr,ki);
}
