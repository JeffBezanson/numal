#include <math.h>
#include <stdio.h>

real_t p(real_t x)
{
	return 1.0;
}

real_t q(real_t x)
{
	return cos(x);
}

real_t r(real_t x)
{
	return exp(x);
}

real_t f(real_t x)
{
	return sin(x)*(1.0+exp(x)+2.0*cos(x));
}

void main ()
{
	void femhermsym(real_t [], real_t [], int, real_t (*)(real_t),
			real_t (*)(real_t), real_t (*)(real_t), real_t (*)(real_t),
			int, real_t []);
	int n,i,order;
	real_t pi,x[11],y[19],e[5],rho1,rho2,d1,d2;

	printf("FEMHERMSYM delivers:\n");
	for (n=5; n<=10; n += 5) {
		e[1]=e[3]=0.0;  e[2]=1.0;  e[4] = -1.0;
		pi=3.14159265358979;
		for (i=0; i<=n; i++) x[i]=pi*i/n;
		printf("N=%2d\n",n);
		for (order=4; order<=8; order += 2) {
			femhermsym(x,y,n,p,q,r,f,order,e);
			rho1=rho2=0.0;
			for (i=1; i<=n-1; i++) {
				d1=fabs(y[2*i-1]-sin(x[i]));
				if (rho1 < d1) rho1=d1;
				d2=fabs(y[2*i]-cos(x[i]));
				if (rho2 < d2) rho2=d2;
			}
			printf("     ORDER=%1d\n"
				"            MAX. ABS(Y[2*I-1]-Y(X[I])) = %7.3e\n"
				"            MAX. ABS(Y[2*I]-Y'(X[I]))  = %7.3e\n",
				order,rho1,rho2);
		}
	}
}

