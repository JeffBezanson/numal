#include "../real.h"


void bessj(real_t x, int n, real_t j[])
{
	if (x == 0.0) {
		j[0]=1.0;
		for (; n>=1; n--) j[n]=0.0;
	} else {
		int start(real_t, int, int);
		int l,m,nu,signx;
		real_t x2,r,s;
		signx = (x > 0.0) ? 1 : -1;
		x=fabs(x);
		r=s=0.0;
		x2=2.0/x;
		l=0;
		nu=start(x,n,0);
		for (m=nu; m>=1; m--) {
			r=1.0/(x2*m-r);
			l=2-l;
			s=r*(l+s);
			if (m <= n) j[m]=r;
		}
		j[0]=r=1.0/(1.0+s);
		for (m=1; m<=n; m++) r = j[m] *= r;
		if (signx < 0.0)
			for (m=1; m<=n; m += 2) j[m] = -j[m];
	}
}
