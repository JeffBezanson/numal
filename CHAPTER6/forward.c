#include "../real.h"
void forward(real_t x, real_t p, real_t q, real_t i0, real_t i1,
				int nmax, real_t i[])
{
	int m,n;
	real_t y,r,s;

	i[0]=i0;
	if (nmax > 0) i[1]=i1;
	m=nmax-1;
	r=p+q-1.0;
	y=1.0-x;
	for (n=1; n<=m; n++) {
		s=(n+r)*y;
		i[n+1]=((n+q+s)*i[n]-s*i[n-1])/(n+q);
	}
}
