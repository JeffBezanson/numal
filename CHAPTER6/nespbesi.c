#include "../real.h"


void nonexpspherbessi(real_t x, int n, real_t i[])
{
	if (x == 0.0) {
		i[0]=1.0;
		for (; n>=1; n--) i[n]=0.0;
	} else {
		int start(real_t, int, int);
		int m;
		real_t x2,r;
		x2=x+x;
		i[0] = x2 = ((x == 0.0) ? 1.0 : ((x2 < 0.7) ?
						sinh(x)/(x*exp(x)) : (1.0-exp(-x2))/x2));
		if (n != 0) {
			r=0.0;
			m=start(x,n,1);
			for (; m>=1; m--) {
				r=1.0/((m+m+1)/x+r);
				if (m <= n) i[m]=r;
			}
			for (m=1; m<=n; m++) x2 = i[m] *= x2;
		}
	}
}
