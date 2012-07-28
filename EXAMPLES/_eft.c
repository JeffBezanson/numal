#include <math.h>
#include <stdio.h>

real_t expt,lnt,u0,u1,u2,accuracy;

void derivative(real_t t, int m0, int m, int i, real_t u[])
{
	if (i == 1) {
		expt=exp(t);
		lnt=log(t);
		u0=u[0];
		u1=u[0]=expt*(lnt-u0)+1.0/t;
	} else if (i == 2)
		u2=u[0]=expt*(lnt-u0-u1+1.0/t)-1.0/t/t;
	else
		u[0]=expt*(lnt-u0-2.0*u1-u2+2.0/t-1.0/t/t)+2.0/t/t/t;
}

real_t sigma(real_t t, int m0, int m)
{
	return exp(t);
}

real_t diameter(real_t t, int m0, int m)
{
	return 2.0*exp(2.0*t/3.0);
}

real_t aeta1(real_t t, int m0, int m)
{
	return accuracy/10.0;
}

real_t aeta2(real_t t, int m0, int m)
{
	return accuracy/10.0*((t < 3.0) ? 1.0 : exp(2.0*(t-3.0)));
}

real_t reta1(real_t t, int m0, int m)
{
	return accuracy;
}

real_t reta2(real_t t, int m0, int m)
{
	return accuracy*((t < 3.0) ? 1.0 : exp(2.0*(t-3.0)));
}

void out(real_t t, real_t te, int m0, int m, real_t u[],
			int k, real_t eta, real_t rho)
{
	if (t == te) printf(" %3d  %e  ",k,u[0]);
}

void main ()
{
	void eft(real_t *, real_t, int, int, real_t [],
			real_t (*)(real_t, int, int), real_t,
			real_t (*)(real_t, int, int),
			void (*)(real_t, int, int, int, real_t []),
			int *, real_t, int,
			real_t (*)(real_t, int, int), real_t (*)(real_t, int, int),
			real_t *, real_t *, real_t, real_t *,
			void (*)(real_t, real_t, int, int, real_t [],
							int, real_t, real_t));
	int j,k,l;
	real_t t,te,te1,te2,eta,rho,pi,hs,u[1];

	printf("The results with EFT are:\n"
		"   K     U(TE1)       K     U(TE2)      RETA\n");
	pi=4.0*atan(1.0);
	te1=exp(1.0);
	te2=exp(2.0);
	accuracy=1.0;
	for (j=1; j<=4; j++) {
		accuracy *= 1.0e-1;
		t=0.01;
		u[0]=log(t);
		k=0;
		hs=0.0;
		for (l=1; l<=2; l++) {
			te = (l == 1) ? te1 : te2;
			eft(&t,te,0,0,u,sigma,pi,diameter,derivative,&k,
				1.5,2,aeta1,reta1,&eta,&rho,1.0e-4,&hs,out);
		}
		printf(" %6.1e\n",accuracy);
	}
	printf("\nWith relaxed accuracy conditions for t > 3 :\n"
		"   K     U(TE1)       K     U(TE2)      RETA\n");
	accuracy=1.0;
	for (j=1; j<=4; j++) {
		accuracy *= 1.0e-1;
		t=0.01;
		u[0]=log(t);
		k=0;
		hs=0.0;
		for (l=1; l<=2; l++) {
			te = (l == 1) ? te1 : te2;
			eft(&t,te,0,0,u,sigma,pi,diameter,derivative,&k,
				1.5,2,aeta2,reta2,&eta,&rho,1.0e-4,&hs,out);
		}
		printf(" %6.1e\n",accuracy);
	}
}

