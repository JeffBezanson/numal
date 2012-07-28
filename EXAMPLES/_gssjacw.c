#include <math.h>
#include <stdio.h>
void main ()
{
	void gssjacwghts(int, real_t, real_t, real_t [], real_t []);
	int n;
	real_t alfa,beta,ind,x[6],w[6];

	alfa=1.0;
	beta=2.0;
	n=5;
	ind=0.0;
	gssjacwghts(n,alfa,beta,x,w);
	for (n=1; n<=5; n++) ind += w[n]*exp(x[n]);
	printf("Delivers:  %13.6e\n",ind-2.0*exp(1.0)+10.0/exp(1.0));
}

