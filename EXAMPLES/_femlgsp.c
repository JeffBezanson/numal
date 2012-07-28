#include <math.h>
#include <stdio.h>

int nc;

real_t r(real_t x)
{
	return 1.0;
}

real_t f(real_t x)
{
	return (12+4*nc)*x*x+1.0-x*x*x*x;
}

void main ()
{
	void femlagspher(real_t [], real_t [], int, int, real_t (*)(real_t),
						real_t (*)(real_t), int, real_t []);
	int n,i,order;
	real_t x[21],y[21],e[7],rho,d;

	printf("FEMLAGSPHER delivers:\n");
	for (n=10; n<=20; n += 10)
		for (nc=0; nc<=2; nc++) {
			e[2]=e[4]=1.0;
			e[1]=e[3]=e[5]=e[6]=0.0;
			for (i=0; i<=n; i++) x[i]=(real_t)(i)/(real_t)(n);
			printf("N= %2d    NC=%2d\n",n,nc);
			for (order=2; order<=4; order += 2) {
				femlagspher(x,y,n,nc,r,f,order,e);
				rho=0.0;
				for (i=0; i<=n; i++) {
					d=fabs(y[i]-1.0+x[i]*x[i]*x[i]*x[i]);
					if (rho < d) rho=d;
				}
				printf("     ORDER=%1d    MAX.ERROR= %7.3e\n",order,rho);
			}
		}
}

