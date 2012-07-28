#include <math.h>
#include <stdio.h>

real_t expt,lnt,c0,c1,c2,c3;

void der(real_t t, int m0, int m, int i, real_t a[])
{
	if (i == 1) {
		expt=exp(t);
		lnt=log(t);
		c0=a[0];
		c1 = a[0] = -expt*c0+1.0/t+expt*lnt;
	}
	if (i == 2) c2=a[0]=expt*(lnt+1.0/t-c0-c1)-1.0/t/t;
	if (i == 3)
		c3=a[0]=expt*(lnt+2.0/t-c0-2.0*c1-c2-1.0/t/t)+2.0/t/t/t;
	if (i == 4)
		a[0]=c3-2.0*(1.0+3.0/t)/t/t/t+expt*((1.0-(2.0-2.0/t)/t)/t-
				c1-c2*2.0-c3);
}

real_t sigma(real_t t, int m0, int m)
{
	return exp(t);
}

real_t aeta(real_t t, int m0, int m)
{
	return 1.0e-5;
}

real_t reta(real_t t, int m0, int m)
{
	return 1.0e-4;
}

void op(real_t t, real_t te, int m0, int m, real_t u[],
			int k, real_t eta, real_t rho)
{
	if (t == te)
		printf("\nNumber of steps: %3d\n"
			"Solution:  T = %8.6f    U(T) = %8.6f\n",k,t,u[0]);
}

void main ()
{
	real_t *allocate_real_vector(int, int);
	void free_real_vector(real_t *, int);
	void modifiedtaylor(real_t *, real_t, int, int, real_t [],
			real_t (*)(real_t, int, int), real_t,
			void (*)(real_t, int, int, int, real_t []),
			int *, real_t [], real_t, int, real_t (*)(real_t, int, int),
			real_t (*)(real_t, int, int), real_t *, real_t *,
			void (*)(real_t, real_t, int, int, real_t [],
						int, real_t, real_t));
	int j,k;
	real_t t,te,eta,rho,u[1],*data;

	data=allocate_real_vector(-2,4);
	printf("The results with MODIFIEDTAYLOR are:\n");
	data[-2]=4.0;  data[-1]=3.0;  data[0]=6.025;  data[1]=1.0;
	data[2]=0.5;  data[3]=1.0/6.0;  data[4]=0.018455702;
	t=u[0]=1.0e-2;
	k=0;
	for (j=1; j<=2; j++) {
		te = (j == 1) ? exp(1.0) : te*te;
		modifiedtaylor(&t,te,0,0,u,sigma,1.0e-4,der,&k,data,1.5,
				1,aeta,reta,&eta,&rho,op);
	}
	free_real_vector(data,-2);
}

