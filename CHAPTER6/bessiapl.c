#include "../real.h"


void bessiaplusn(real_t a, real_t x, int n, real_t ia[])
{
	if (x == 0.0) {
		ia[0] = (a == 0.0) ? 1.0 : 0.0;
		for (; n>=1; n--) ia[n]=0.0;
	} else if (a == 0.0) {
		void bessi(real_t, int, real_t []);
		bessi(x,n,ia);
	} else if (a == 0.5) {
		void nonexpspherbessi(real_t, int, real_t []);
		real_t c;
		c=0.797884560802865*sqrt(fabs(x))*exp(fabs(x));
		nonexpspherbessi(x,n,ia);
		for (; n>=0; n--) ia[n] *= c;
	} else {
		void nonexpbessiaplusn(real_t, real_t, int, real_t[]);
		real_t expx;
		expx=exp(fabs(x));
		nonexpbessiaplusn(a,x,n,ia);
		for (; n>=0; n--) ia[n] *= expx;
	}
}
