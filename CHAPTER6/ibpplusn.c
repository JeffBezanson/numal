#include "../real.h"
void ibpplusn(real_t x, real_t p, real_t q, int nmax, real_t eps,
					real_t i[])
{
	void ixqfix(real_t, real_t, real_t, int, real_t, real_t []);
	void ixpfix(real_t, real_t, real_t, int, real_t, real_t []);
	int n;

	if (x == 0.0 || x == 1.0)
		for (n=0; n<=nmax; n++) i[n]=x;
	else {
		if (x <= 0.5)
			ixqfix(x,p,q,nmax,eps,i);
		else {
			ixpfix(1.0-x,q,p,nmax,eps,i);
			for (n=0; n<=nmax; n++) i[n]=1.0-i[n];
		}
	}
}
