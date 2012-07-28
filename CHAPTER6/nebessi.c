#include "../real.h"


void nonexpbessi(real_t x, int n, real_t i[])
{
	if (x == 0.0) {
		i[0]=1.0;
		for (; n>=1; n--) i[n]=0.0;
	} else {
		int start(real_t, int, int);
		int k,negative;
		real_t x2,r,s;
		negative = (x < 0.0);
		x=fabs(x);
		r=s=0.0;
		x2=2.0/x;
		k=start(x,n,1);
		for (; k>=1; k--) {
			r=1.0/(r+x2*k);
			s=r*(2.0+s);
			if (k <= n) i[k]=r;
		}
		i[0]=r=1.0/(1.0+s);
		if (negative)
			for (k=1; k<=n; k++) r = i[k] *= (-r);
		else
			for (k=1; k<=n; k++) r = i[k] *= r;
	}
}
