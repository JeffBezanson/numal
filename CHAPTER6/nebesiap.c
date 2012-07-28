#include "../real.h"


void nonexpbessiaplusn(real_t a, real_t x, int n, real_t ia[])
{
	if (x == 0.0) {
		ia[0] = (a == 0.0) ? 1.0 : 0.0;
		for (; n>=1; n--) ia[n]=0.0;
	} else if (a == 0.0) {
		void nonexpbessi(real_t, int, real_t []);
		nonexpbessi(x,n,ia);
	} else if (a == 0.5) {
		void nonexpspherbessi(real_t, int, real_t []);
		real_t c;
		c=0.797884560802865*sqrt(x);
		nonexpspherbessi(x,n,ia);
		for (; n>=0; n--) ia[n] *= c;
	} else {
		real_t numal_gamma(real_t);
		int start(real_t, int, int);
		int m,nu;
		real_t r,s,labda,l,a2,x2;
		a2=a+a;
		x2=2.0/x;
		l=1.0;
		nu=start(x,n,1);
		r=s=0.0;
		for (m=1; m<=nu; m++) l=l*(m+a2)/(m+1);
		for (m=nu; m>=1; m--) {
			r=1.0/(x2*(a+m)+r);
			l=l*(m+1)/(m+a2);
			labda=l*(m+a)*2.0;
			s=r*(labda+s);
			if (m <= n) ia[m]=r;
		}
		ia[0]=r=1.0/(1.0+s)/numal_gamma(1.0+a)/pow(x2,a);
		for (m=1; m<=n; m++) r = ia[m] *= r;
	}
}
