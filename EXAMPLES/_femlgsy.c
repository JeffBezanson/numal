#include <math.h>
#include <stdio.h>

real_t r(real_t x)
{
	return cos(x);
}

real_t p(real_t x)
{
	return exp(x);
}

real_t f(real_t x)
{
	return exp(x)*(sin(x)-cos(x))+sin(2.0*x)/2.0;
}

void main ()
{
	void femlagsym(real_t [], real_t [], int, real_t (*)(real_t),
					real_t (*)(real_t), real_t (*)(real_t), int, real_t []);
	int n,i,order;
	real_t pi,x[21],y[21],e[7],rho,d;

	printf("FEMLAGSYM delivers:\n");
	for (n=10; n<=20; n += 10) {
		e[1]=e[4]=1.0;
		e[2]=e[3]=e[5]=e[6]=0.0;
		pi=3.14159265358979;
		for (i=0; i<=n; i++) x[i]=pi*i/n;
		printf("N=%2d\n",n);
		for (order=2; order<=6; order += 2) {
			femlagsym(x,y,n,p,r,f,order,e);
			rho=0.0;
			for (i=0; i<=n; i++) {
				d=fabs(y[i]-sin(x[i]));
				if (rho < d) rho=d;
			}
			printf("     ORDER=%1d    MAX.ERROR= %7.3e\n",order,rho);
		}
	}
}

