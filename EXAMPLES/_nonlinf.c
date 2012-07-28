#include <math.h>
#include <stdio.h>

int nc;

real_t f(real_t x, real_t y, real_t z)
{
	return exp(y)+exp(z)-exp(1.0-x*x)-exp(-2.0*x)-2.0-2*nc;
}

real_t fy(real_t x, real_t y, real_t z)
{
	return exp(y);
}

real_t fz(real_t x, real_t y, real_t z)
{
	return exp(z);
}

void main ()
{
	void nonlinfemlagskew(real_t [], real_t [], int,
				real_t (*)(real_t, real_t, real_t),
				real_t (*)(real_t, real_t, real_t),
				real_t (*)(real_t, real_t, real_t), int, real_t []);
	int n,i;
	real_t x[51],y[51],e[7],rho,d;

	printf("NONLINFEMLAGSKEW delivers:\n");
	for (nc=0; nc<=2; nc++)
		for (n=25; n<=50; n += 25) {
			e[2]=e[4]=1.0;
			e[1]=e[3]=e[5]=e[6]=0.0;
			for (i=0; i<=n; i++) {
				x[i]=(real_t)(i)/(real_t)(n);
				y[i]=0.0;
			}
			printf(" N=%2d    NC=%1d",n,nc);
			nonlinfemlagskew(x,y,n,f,fy,fz,nc,e);
			rho=0.0;
			for (i=0; i<=n; i++) {
				d=fabs(y[i]-1.0+x[i]*x[i]);
				if (rho < d) rho=d;
				}
			printf("    MAX.ERROR= %7.3e\n",rho);
		}
}

